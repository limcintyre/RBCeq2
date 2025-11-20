#!/usr/bin/env python3

import argparse
import sys
from collections import defaultdict
from functools import partial
from multiprocessing import Pool
from pathlib import Path
from typing import Callable
import os
import pandas as pd
from icecream import ic
from loguru import logger
from typing import Mapping

import rbceq2.core_logic.co_existing as co
import rbceq2.core_logic.data_procesing as dp
import rbceq2.filters.geno as filt
import rbceq2.filters.phased as filt_phase
import rbceq2.filters.knops as filt_co
import rbceq2.phenotype.choose_pheno as ph
from rbceq2.core_logic.constants import PhenoType, DB_VERSION, VERSION
from rbceq2.core_logic.utils import compose, get_allele_relationships
from rbceq2.db.db import (
    Db,
    prepare_db,
    DbDataConsistencyChecker,
    build_antigen_map_for_checks,
)
from rbceq2.IO.PDF_reports import generate_all_reports
from rbceq2.IO.record_data import (
    check_VCF,
    configure_logging,
    log_validation,
    record_filtered_data,
    save_df,
    stamps,
)
from rbceq2.IO.vcf import (
    VCF,
    filter_VCF_to_BG_variants,
    read_vcf,
    check_if_multi_sample_vcf,
    split_vcf_to_dfs,
    build_intervals,
)

from rbceq2.core_logic.large_variants import (
    SnifflesVcfSvReader,
    load_db_defs,
    SvMatcher,
    select_best_per_vcf,
)


def parse_args(args: list[str]) -> argparse.Namespace:
    """Parse command-line arguments for somatic variant calling.

    Args:
        args (List[str]): List of strings representing the command-line arguments.

    Returns:
        argparse.Namespace: An object containing the parsed command-line options.

    This function configures and interprets command-line options for a somatic
    variant caller. It expects paths to VCF files, a database file, and allows
    specification of output options and genomic references.
    """
    parser = argparse.ArgumentParser(
        description="Calls ISBT defined alleles from VCF/s",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        usage="rbceq2 --vcf example.vcf.gz --out example --reference_genome GRCh37",
    )
    version_str = f"%(prog)s {VERSION} (DB: {DB_VERSION})"

    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=version_str,
        help="Show program's version number and exit.",
    )
    parser.add_argument(
        "--vcf",
        type=lambda p: Path(p).absolute(),
        help=(
            "Path to VCF file/s. Give a folder if you want to pass multiple "
            "separate files (file names must end in .vcf or .vcf.gz), or "
            "alternatively give a file if using a multi-sample VCF."
        ),
    )
    parser.add_argument(
        "--out", type=lambda p: Path(p).absolute(), help="Prefix for output files"
    )
    parser.add_argument(
        "--no_filter",
        action="store_true",
        help="Use all variants, not just those where FILTER = PASS in the VCF",
        default=False,
    )
    parser.add_argument(
        "--processes",
        type=int,
        help=(
            "Number of processes. I.e., how many CPUs are available? ~1GB RAM required "
            "per process"
        ),
        default=1,
    )
    parser.add_argument(
        "--reference_genome",
        type=str,
        help=("GRCh37/8"),
        choices=["GRCh37", "GRCh38"],
        required=True,
    )
    parser.add_argument(
        "--phased",
        action="store_true",
        help="Use phase information",
        default=False,
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug logging. If not set, logging will be at info level.",
        default=False,
    )
    parser.add_argument(
        "--validate",
        action="store_true",
        help="Enable VCF validation. Doubles run time. Might help you identify input issues",
        default=False,
    )
    parser.add_argument(
        "--PDFs",
        action="store_true",
        help="Generate a per sample PDF report",
        default=False,
    )
    parser.add_argument(
        "--HPAs",
        action="store_true",
        help="Generate results for HPA",
        default=False,
    )
    parser.add_argument(
        "--min_size",
        type=int,
        help=("Minimum size indel/SV to apply fuzzy matching to"),
        default=10,
    )
    parser.add_argument(
        "--RH",
        action="store_true",
        help="Generate results for RHD and RHCE. WARNING! Long read only!",
        default=False,
    )

    return parser.parse_args(args)


def main():
    ic("Running RBCeq2...")

    start = pd.Timestamp.now()
    args = parse_args(sys.argv[1:])
    exclude = ["C4A", "C4B", "CD99", "ATP11C"]  # , "RHD", "RHCE"]
    if not args.RH:
        exclude += ["RHD", "RHCE"]
    if not args.HPAs:
        exclude += [f"HPA{i}" for i in range(50)]
    # Configure logging
    UUID = configure_logging(args)

    logger.debug("Logger configured for debug mode.")
    logger.info("Application started.")

    # 1. Prepare the DataFrame
    logger.info("Preparing database DataFrame...")
    db_df = prepare_db()
    logger.info("Database DataFrame prepared.")

    # 2. Run consistency checks on the prepared DataFrame
    DbDataConsistencyChecker.run_all_checks(df=db_df)
    # If any check fails, an exception will be raised here, and the program will halt.

    # 3. If all checks pass, proceed to create the Db object
    logger.info("Consistency checks passed. Initializing Db object...")
    db = Db(ref=args.reference_genome, df=db_df)
    logger.info("Db object initialized.")

    if args.vcf.is_dir():
        patterns = ["*.vcf", "*.vcf.gz"]
        vcfs = [file for pattern in patterns for file in args.vcf.glob(pattern)]
        logger.info(f"{len(vcfs)} single sample VCF/s passed")
        if args.validate:
            with Pool(processes=int(args.processes)) as pool:
                for result_valid, file_name in pool.imap_unordered(
                    check_VCF, list(vcfs)
                ):
                    log_validation(result_valid, file_name)
    else:
        if args.validate:
            result_valid, file_name = check_VCF(args.vcf)
            log_validation(result_valid, file_name)
        actually_multi_vcf = check_if_multi_sample_vcf(args.vcf)
        if actually_multi_vcf:
            intervals = build_intervals(db_df, args.reference_genome)
            multi_vcf = read_vcf(str(args.vcf), intervals)
            logger.info("Multi sample VCF passed")
            filtered_multi_vcf = filter_VCF_to_BG_variants(
                multi_vcf, db.unique_variants
            )
            vcfs = split_vcf_to_dfs(filtered_multi_vcf)
            time_str = stamps(start)
            logger.info(f"VCFs loaded in {time_str}")
            print(f"VCFs loaded in {time_str}")
        else:
            logger.info("1 single sample VCF passed")
            vcfs = [args.vcf]
    mapping = build_antigen_map_for_checks(db_df)

    all_alleles = defaultdict(list)
    for a in db.make_alleles():
        all_alleles[a.blood_group].append(a)
    allele_relationships = get_allele_relationships(all_alleles, int(args.processes))
    dfs_geno = {}
    dfs_pheno_numeric = {}
    dfs_pheno_alphanumeric = {}
    with Pool(processes=int(args.processes)) as pool:
        find_hits_db = partial(
            find_hits,
            db,
            args=args,
            allele_relationships=allele_relationships,
            excluded=exclude,
            ant_mapping=mapping,
        )
        for results in pool.imap_unordered(find_hits_db, list(vcfs)):
            if results is not None:
                sample, genos, numeric_phenos, alphanumeric_phenos, _, _ = results
                dfs_geno[sample] = genos
                dfs_pheno_numeric[sample] = numeric_phenos
                dfs_pheno_alphanumeric[sample] = alphanumeric_phenos
                record_filtered_data(results)
                sep = "##############"
                logger.debug(f"\n {sep} End log for sample: {sample} {sep}\n")

    df_geno = pd.DataFrame.from_dict(dfs_geno, orient="index")
    df_geno = df_geno.replace("", "Undetermined/Undetermined")
    save_df(df_geno, f"{args.out}_geno.tsv", UUID)
    df_pheno_numeric = pd.DataFrame.from_dict(dfs_pheno_numeric, orient="index")
    save_df(df_pheno_numeric, f"{args.out}_pheno_numeric.tsv", UUID)
    df_pheno_alpha = pd.DataFrame.from_dict(dfs_pheno_alphanumeric, orient="index")
    save_df(df_pheno_alpha, f"{args.out}_pheno_alphanumeric.tsv", UUID)
    if args.PDFs:
        generate_all_reports(df_geno, df_pheno_alpha, df_pheno_numeric, args.out, UUID)

    time_str = stamps(start)
    logger.info(f"{len(dfs_geno)} VCFs processed in {time_str}")
    print(
        f"\n✅ Complete! {len(dfs_geno)} VCFs processed in {time_str}. \n💾 Results saved successfully."
    )


def find_hits(
    db: Db,
    vcf: tuple[pd.DataFrame, str] | Path,
    args: argparse.Namespace,
    allele_relationships: dict[str, dict[str, bool]],
    excluded: list[str],
    ant_mapping: Mapping[str, Mapping[str, str]],
) -> pd.DataFrame | None:
    intervals = build_intervals(db.df, args.reference_genome)
    if isinstance(vcf, Path):
        vcf_filtered_by_500kb_padded_bed = read_vcf(str(vcf), intervals)
        vcf = VCF(
            vcf_filtered_by_500kb_padded_bed,
            db.lane_variants,
            db.unique_variants,
            vcf.stem,
        )
    else:
        vcf = VCF(vcf, db.lane_variants, db.unique_variants, vcf[-1])
    reader = SnifflesVcfSvReader(df=vcf.df, min_size=args.min_size)
    events = list(reader.events())

    db_defs = load_db_defs(db.df)
    matcher = SvMatcher()
    matches = matcher.match(db_defs, events)
    best = select_best_per_vcf(matches, tie_tol=1e-9)
    var_map = {}

    if best:
        for match in best:
            vcf.variants[f"{match.vcf.chrom}:{match.db.raw}"] = dict(
                zip(match.vcf.sample_fmt.split(":"), match.vcf.sample_value.split(":"))
            )
            var_map[f"{match.vcf.chrom}:{match.db.raw}"] = match.variant
    res = dp.raw_results(db, vcf, excluded, var_map, matches)
    res = dp.make_blood_groups(res, vcf.sample)

    pipe: list[Callable] = [
        partial(
            dp.only_keep_alleles_if_FILTER_PASS, df=vcf.df, no_filter=args.no_filter
        ),
        partial(dp.make_variant_pool, vcf=vcf),
        partial(dp.modify_variant_pool_if_large_indel),
        partial(dp.modify_allele_pool_if_large_indel),
        partial(
            dp.add_phasing,
            phased=args.phased,
            variant_metrics=vcf.variants,
            phase_sets=vcf.phase_sets,
        ),
        partial(dp.modify_phase_of_large_indel, phased=args.phased),
        partial(dp.modify_variant_phase_pool_if_large_indel),
        partial(filt_phase.remove_unphased, phased=args.phased),
        partial(dp.process_genetic_data, reference_alleles=db.reference_alleles),
        partial(
            dp.find_what_was_excluded_due_to_rank,
            reference_alleles=db.reference_alleles,
        ),
        filt.cant_not_include_null,
        partial(
            filt.filter_pairs_on_antithetical_zygosity, antitheticals=db.antitheticals
        ),
        partial(
            filt_phase.filter_pairs_by_phase,
            phased=args.phased,
            reference_alleles=db.reference_alleles,
        ),
        partial(
            filt_phase.no_defining_variant,
            phased=args.phased,
        ),
        partial(
            filt_phase.ref_not_phased,
            phased=args.phased,
        ),
        partial(
            filt_phase.cant_be_hom_ref_due_to_HET_SNP,
            phased=args.phased,
        ),
        co.homs,
        co.max_rank,
        partial(co.prep_co_putative_combos, allele_relationships=allele_relationships),
        co.add_co_existing_alleles,
        partial(
            co.add_co_existing_allele_and_ref, reference_alleles=db.reference_alleles
        ),
        co.filter_redundant_pairs,
        co.mush,
        partial(
            co.list_excluded_co_existing_pairs, reference_alleles=db.reference_alleles
        ),
        partial(
            filt_co.filter_coexisting_pairs_on_antithetical_zygosity,
            antitheticals=db.antitheticals,
        ),
        partial(filt_co.remove_unphased_co, phased=args.phased),
        filt.cant_pair_with_ref_cuz_trumped,
        partial(
            filt.antithetical_modifying_SNP_is_HOM,
            antitheticals=db.antitheticals,
        ),
        filt.ensure_HET_SNP_used,
        filt.ABO_cant_pair_with_ref_cuz_261delG_HET,
        filt.cant_pair_with_ref_cuz_SNPs_must_be_on_other_side,
        filt.filter_HET_pairs_by_weight,
        filt.filter_pairs_by_context,
        filt.impossible_alleles,
        partial(filt_phase.impossible_alleles_phased, phased=args.phased),
        partial(
            filt_phase.filter_if_all_HET_vars_on_same_side_and_phased,
            phased=args.phased,
        ),
        partial(
            filt_phase.filter_on_in_relationship_if_HET_vars_on_dif_side_and_phased,
            phased=args.phased,
        ),
        partial(
            filt_phase.filter_on_in_relationship_if_all_HOM_and_phased,
            phased=args.phased,
        ),
        partial(
            filt_phase.filter_on_in_relationship_when_HOM_cant_be_on_one_side,
            phased=args.phased,
        ),
        partial(filt_phase.rm_ref_if_2x_HET_phased, phased=args.phased),
        partial(filt_phase.low_weight_hom, phased=args.phased),
        filt_co.ensure_co_existing_HET_SNP_used,
        filt_co.filter_co_existing_pairs,
        filt_co.filter_co_existing_in_other_allele,
        filt_co.filter_co_existing_with_normal,  # has to be after normal filters!!!!!!!
        filt_co.filter_co_existing_subsets,
        dp.get_genotypes,
        dp.add_CD_to_XG,
    ]
    preprocessor = compose(*pipe)
    res = preprocessor(res)

    res = dp.add_refs(db, res, excluded)

    # merge FUT 1 and 2
    fut2s = res["FUT2"].genotypes.copy()
    fut1s = res["FUT1"].genotypes.copy()
    for allele_pair in fut2s:
        res["FUT1"].genotypes.append(allele_pair)
    for allele_pair in fut1s:
        res["FUT2"].genotypes.append(allele_pair)

    formated_called_genos = {k: ",".join(bg.genotypes) for k, bg in res.items()}
    pipe2: list[Callable] = [
        partial(ph.add_ref_phenos, df=db.df),
        partial(ph.instantiate_antigens, ant_type=PhenoType.numeric),
        partial(ph.instantiate_antigens, ant_type=PhenoType.alphanumeric),
        partial(ph.get_phenotypes1, ant_type=PhenoType.numeric),
        partial(ph.get_phenotypes1, ant_type=PhenoType.alphanumeric),
        partial(ph.get_phenotypes2, ant_type=PhenoType.numeric),
        partial(ph.get_phenotypes2, ant_type=PhenoType.alphanumeric),
        partial(ph.internal_anithetical_consistency_HET, ant_type=PhenoType.numeric),
        partial(
            ph.internal_anithetical_consistency_HET, ant_type=PhenoType.alphanumeric
        ),
        partial(ph.internal_anithetical_consistency_HOM, ant_type=PhenoType.numeric),
        partial(
            ph.internal_anithetical_consistency_HOM, ant_type=PhenoType.alphanumeric
        ),
        partial(ph.include_first_antithetical_pair, ant_type=PhenoType.numeric),
        partial(ph.include_first_antithetical_pair, ant_type=PhenoType.alphanumeric),
        partial(ph.sort_antigens, ant_type=PhenoType.numeric),
        partial(ph.sort_antigens, ant_type=PhenoType.alphanumeric),
        partial(ph.phenos_to_str, ant_type=PhenoType.numeric),
        partial(ph.phenos_to_str, ant_type=PhenoType.alphanumeric),
        partial(ph.modify_FY, ant_type=PhenoType.numeric),
        partial(ph.compare_numeric_ants_to_alphanumeric, mapping=ant_mapping),
        ph.combine_anitheticals,
        partial(ph.modify_FY, ant_type=PhenoType.alphanumeric),
        partial(ph.modify_KEL, ant_type=PhenoType.alphanumeric),
        partial(ph.modify_CROM, ant_type=PhenoType.alphanumeric),
        partial(ph.re_order_KEL, ant_type=PhenoType.alphanumeric),
        partial(ph.modify_MNS, ant_type=PhenoType.alphanumeric),
        partial(ph.modify_FY2, ant_type=PhenoType.alphanumeric),
        partial(ph.modify_RHD, ant_type=PhenoType.numeric),
        partial(ph.modify_RHD, ant_type=PhenoType.alphanumeric),
        
    ]

    preprocessor2 = compose(*pipe2)
    res = preprocessor2(res)
    res = ph.FUT3(res)
    res = ph.FUT1(res)

    formated_called_numeric_phenos = {
        k: " | ".join(sorted(set(bg.phenotypes[PhenoType.numeric].values())))
        for k, bg in res.items()
    }
    formated_called_alphanumeric_phenos = {
        k: " | ".join(sorted(set(bg.phenotypes[PhenoType.alphanumeric].values())))
        for k, bg in res.items()
    }

    return (
        vcf.sample,
        formated_called_genos,
        formated_called_numeric_phenos,
        formated_called_alphanumeric_phenos,
        res,
        var_map,
    )


if __name__ == "__main__":
    main()
