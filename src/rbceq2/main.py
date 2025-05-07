#!/usr/bin/env python3

import argparse
import sys
from collections import defaultdict
from functools import partial
from multiprocessing import Pool
from pathlib import Path
from typing import Callable

import pandas as pd
from icecream import ic
from loguru import logger

import rbceq2.core_logic.co_existing as co
import rbceq2.core_logic.data_procesing as dp
import rbceq2.filters.geno as gf
import rbceq2.phenotype.choose_pheno as ph
from rbceq2.core_logic.constants import PhenoType
from rbceq2.core_logic.utils import compose, get_allele_relationships
from rbceq2.db.db import Db, prepare_db, DbInternalConsistencyCheck
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
        "--depth", type=int, help="Minimum number of reads for a variant", default=10
    )
    parser.add_argument(
        "--quality",
        type=int,
        help="Minimum average genotype quality for a variant",
        default=10,
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
        "--microarray",
        action="store_true",
        help="Input is from a microarray.",
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
        help="Enable VCF validation. Doubles run time",
        default=False,
    )
    parser.add_argument(
        "--PDFs",
        action="store_true",
        help="Generate a per sample PDF report",
        default=False,
    )

    return parser.parse_args(args)


def main():
    ic("Running RBCeq2...")

    start = pd.Timestamp.now()
    args = parse_args(sys.argv[1:])

    # Configure logging
    UUID = configure_logging(args)

    logger.debug("Logger configured for debug mode.")
    logger.info("Application started.")
    db_df = prepare_db()
    #db_check = DbInternalConsistencyCheck(ref=args.reference_genome, df = db_df)
    db = Db(ref=args.reference_genome, df = db_df)

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
            multi_vcf = read_vcf(args.vcf)
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
        )
        for results in pool.imap_unordered(find_hits_db, list(vcfs)):
            if results is not None:
                sample, genos, numeric_phenos, alphanumeric_phenos, res = results
                dfs_geno[sample] = genos
                dfs_pheno_numeric[sample] = numeric_phenos
                dfs_pheno_alphanumeric[sample] = alphanumeric_phenos
                record_filtered_data(results)
                sep = "##############"
                logger.debug(f"\n {sep} End log for sample: {sample} {sep}\n")

    df_geno = pd.DataFrame.from_dict(dfs_geno, orient="index")
    save_df(df_geno, f"{args.out}_geno.tsv", UUID)
    df_pheno_numeric = pd.DataFrame.from_dict(dfs_pheno_numeric, orient="index")
    save_df(df_pheno_numeric, f"{args.out}_pheno_numeric.tsv", UUID)
    df_pheno_alpha = pd.DataFrame.from_dict(dfs_pheno_alphanumeric, orient="index")
    save_df(df_pheno_alpha, f"{args.out}_pheno_alphanumeric.tsv", UUID)
    if args.PDFs:
        generate_all_reports(df_geno, df_pheno_alpha, df_pheno_numeric, args.out, UUID)

    time_str = stamps(start)
    logger.info(f"{len(dfs_geno)} VCFs processed in {time_str}")
    print(f"{len(dfs_geno)} VCFs processed in {time_str}")


def find_hits(
    db: Db,
    vcf: tuple[pd.DataFrame, str] | Path,
    args: argparse.Namespace,
    allele_relationships: dict[str, dict[str, bool]],
) -> pd.DataFrame | None:
    vcf = VCF(vcf, db.lane_variants, db.unique_variants)

    res = dp.raw_results(db, vcf)
    res = dp.make_blood_groups(res, vcf.sample)

    pipe: list[Callable] = [
        partial(
            dp.remove_alleles_with_low_read_depth,
            variant_metrics=vcf.variants,
            min_read_depth=args.depth,
            microarray=args.microarray,
        ),
        partial(
            dp.remove_alleles_with_low_base_quality,
            variant_metrics=vcf.variants,
            min_base_quality=args.quality,
            microarray=args.microarray,
        ),
        partial(dp.make_variant_pool, vcf=vcf),
        partial(dp.add_phasing, phased=args.phased, variant_metrics=vcf.variants),
        partial(gf.remove_unphased, phased=args.phased),
        dp.rule_out_impossible_alleles,
        partial(dp.process_genetic_data, reference_alleles=db.reference_alleles),
        partial(
            dp.find_what_was_excluded_due_to_rank,
            reference_alleles=db.reference_alleles,
        ),
        partial(
            gf.filter_pairs_on_antithetical_zygosity, antitheticals=db.antitheticals
        ),
        partial(
            gf.filter_pairs_on_antithetical_modyfying_SNP,
            antitheticals=db.antitheticals,
        ),
        partial(
            gf.filter_pairs_by_phase,
            phased=args.phased,
            reference_alleles=db.reference_alleles,
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
            gf.filter_coexisting_pairs_on_antithetical_zygosity,
            antitheticals=db.antitheticals,
        ),
        gf.cant_pair_with_ref_cuz_trumped,
        gf.ABO_cant_pair_with_ref_cuz_261delG_HET,
        gf.cant_pair_with_ref_cuz_SNPs_must_be_on_other_side,
        gf.filter_HET_pairs_by_weight,
        gf.filter_pairs_by_context,
        gf.ensure_co_existing_HET_SNP_used,
        gf.filter_co_existing_pairs,
        gf.filter_co_existing_in_other_allele,
        gf.filter_co_existing_with_normal,  # has to be after normal filters!!!!!!!
        gf.filter_co_existing_subsets,
        dp.get_genotypes,
        dp.add_CD_to_XG,
    ]
    preprocessor = compose(*pipe)
    res = preprocessor(res)

    res = dp.add_refs(db, res)

    for allele_pair in res["FUT2"].genotypes:
        res["FUT1"].genotypes.append(allele_pair)
    for allele_pair in res["FUT1"].genotypes:
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
        ph.combine_anitheticals,
        partial(ph.modify_FY, ant_type=PhenoType.alphanumeric),
        partial(ph.modify_KEL, ant_type=PhenoType.alphanumeric),
        partial(ph.re_order_KEL, ant_type=PhenoType.alphanumeric),
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
    )


if __name__ == "__main__":
    main()
