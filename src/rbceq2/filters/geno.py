from __future__ import annotations

import operator
from collections import defaultdict
from rbceq2.core_logic.alleles import Allele, BloodGroup, Pair
from rbceq2.core_logic.constants import (
    ABO_DELG_VARIANTS,
    LOW_WEIGHT,
    UNDETERMINED_SLOT,
    AlleleState,
)
from rbceq2.core_logic.utils import (
    BeyondLogicError,
    Zygosity,
    apply_to_dict_values,
    check_available_variants,
    one_HET_variant,
)
from rbceq2.filters.shared_filter_functionality import (
    flatten_alleles,
    proceed,
    all_hom,
    check_var,
)

def split_pair_by_ref(pair: Pair) -> tuple[Allele, Allele]:
    """Split a pair of alleles into reference and non-reference.

    This function assumes that exactly one of the alleles in the pair is marked as a
    reference.

    Args:
        pair (Pair): The pair of alleles to split.

    Returns:
        tuple[Allele, Allele]: A tuple where the first element is the reference allele
                                and the second is the non-reference allele.

    Raises:
        ValueError: If both or neither alleles are marked as reference.
    """

    if pair.allele1.reference and not pair.allele2.reference:
        ref, allele = pair.allele1, pair.allele2
    elif not pair.allele1.reference and pair.allele2.reference:
        ref, allele = pair.allele2, pair.allele1
    else:
        raise ValueError("Both ref")

    return ref, allele


@apply_to_dict_values
def filter_pairs_on_antithetical_zygosity(
    bg: BloodGroup, antitheticals: dict[str, list[str]]
) -> BloodGroup:
    """Process genetic data to identify alleles and genotypes.

    Args:
    ----------
    res (dict[str, list[Allele]]):
        A dictionary mapping genotypes to lists of Allele objects.
    variant_pool_numeric (dict[str, int]):
        A dictionary mapping variants to their counts.

    Returns:
    ----------
    list[str]:
        A list of called genotypes based on processed data.

    Example:
    ----------

    hashed ones not possible because loci that defines the antithetical antigens is HET,
    meaning there will def be 1 FY*01* and 1 FY*02*
    (and in this case FY*02/FY*01 not possible as a modifying SNP exists)
    'FY*01N.01/FY*02'
    #'FY*01N.01/FY*01' not possible
    'FY*02N.01/FY*01'
    #'FY*02/FY*01'
    """

    def find_anti_variant_pair(data_list: list[str]) -> tuple[str] | None:
        """
        Parses a list of comma-separated variant strings to find a matching pair.

        A matching pair consists of two variants that share the same chromosomal
        position but have different alleles (one '_ref' and one alternate).
        These positions act as an antithetical switch between sub types

        Args:
            data_list: A list of strings, where each string can contain one or
                    more comma-separated variants.
                    e.g., ['25408711_ref,25390874_ref', '25408711_G_A']

        Returns:
            A tuple containing the reference variant and the alternate variant
            (e.g., ('25408711_ref', '25408711_G_A')).
            Returns None if no such pair is found.
        """
        variants_by_pos = defaultdict(set)
        for item in data_list:
            variants = item.split(",")
            for variant in variants:
                if not variant:
                    continue
                position = variant.split("_")[0]
                variants_by_pos[position].add(variant)

        for position, unique_variants in variants_by_pos.items():
            if len(unique_variants) == 2:
                ref_variant = None
                alt_variant = None
                for variant in unique_variants:
                    if variant.endswith("_ref"):
                        ref_variant = variant
                    else:
                        alt_variant = variant

                if ref_variant and alt_variant:
                    return (ref_variant, alt_variant)
        return None

    to_remove = []
    if bg.type in antitheticals:
        anti1, anti2 = find_anti_variant_pair(antitheticals[bg.type])
        var_pool_no_chrom = {k.split(":")[1]: v for k, v in bg.variant_pool.items()}
        if var_pool_no_chrom.get(anti1) == Zygosity.HOM:
            return bg
        if var_pool_no_chrom.get(anti2) == Zygosity.HOM:
            return bg
        if anti1 not in var_pool_no_chrom or anti2 not in var_pool_no_chrom:
            # if Lane var doesn't get PASS FILTER
            return bg
        if (
            var_pool_no_chrom.get(anti1) is None
            and var_pool_no_chrom.get(anti2) is None
        ):
            return bg  # KN has alleles without either antithetical SNV

        # remove pairs of same subtype
        for pair in bg.alleles[AlleleState.NORMAL]:
            if pair.allele1.sub_type == pair.allele2.sub_type:
                to_remove.append(pair)
        if to_remove:
            bg.remove_pairs(to_remove, "filter_pairs_on_antithetical_zygosity")

    return bg


@apply_to_dict_values
def antithetical_modifying_SNP_is_HOM(
    bg: BloodGroup, antitheticals: dict[str, list[str]]
) -> BloodGroup:
    """Process genetic data to identify alleles and genotypes.

    Parameters:
    ----------
    res (dict[str, list[Allele]]):
        A dictionary mapping genotypes to lists of Allele objects.
    variant_pool_numeric (dict[str, int]):
        A dictionary mapping variants to their counts.

    Returns
    -------
    list[str]:
        A list of called genotypes based on processed data.

    Example:
    ----------

    'LU*01.19/LU*02' not possible because modifying SNP (45322744_A_G) is hom
    Allele(genotype='LU*02',
            phenotype='LU:2',
            genotype_alt='LU*B',
            phenotype_alt='Lu(a-b+)',
            sample='',
            defining_variants=frozenset({'19:45315445_ref'}),
            weight_geno=LOW_WEIGHT,
            weight_pheno=2,
            reference=True,
            Sub_type='LU*02'),
    Allele(genotype='LU*02.19',
            phenotype='LU:-18,19',
            genotype_alt='.',
            phenotype_alt='Au(a-b+)',
            sample='',
            defining_variants=frozenset({'19:45315445_ref',
                                        '19:45322744_A_G'}), hom
            weight_geno=LOW_WEIGHT,
            weight_pheno=1,
            reference=False,
            Sub_type='LU*02')]},
    sample='128',
    variant_pool={'19:45315445_G_A': 'Heterozygous',
                '19:45315445_ref': 'Heterozygous',
                '19:45322744_A_G': 'Homozygous'},
    genotypes=['LU*01.19/LU*02', 'LU*01.19/LU*02.19'],
    """
    for allele_state in [AlleleState.NORMAL, AlleleState.CO]:
        if not proceed(bg, allele_state):
            continue
        to_remove = []
        if bg.type in antitheticals:
            modifying_SNP = None
            flattened_alleles = flatten_alleles(bg.alleles[allele_state])
            d = defaultdict(set)
            for allele in flattened_alleles:
                if allele.number_of_defining_variants > 1:
                    for variant in allele.defining_variants:
                        if variant.split(":")[1] not in antitheticals[bg.type]:
                            d[allele.sub_type].add(variant)

            if len(d) > 1:
                putative_mod_SNPs = set.union(*d.values())
                if len(putative_mod_SNPs) == 1:
                    modifying_SNP = putative_mod_SNPs.pop()
            if (
                modifying_SNP is not None
                and bg.variant_pool[modifying_SNP] == Zygosity.HOM
            ):
                if len(d) != 2:
                    # Deliberately placed after the homozygosity test rather than before
                    # it. An antithetical system having more than two subtypes is not by
                    # itself a problem - GYPB has three and it is fine, because the single
                    # modifying SNP test below normally rules the case out and nothing
                    # happens. It is only unclear when a group with more than two subtypes
                    # *also* reaches the point of removing pairs, and that has never been
                    # seen. This carries the context because when it does happen the error
                    # is the entire evidence: there is no known case to reproduce from, so
                    # whoever hits it first has to be able to hand it over as it stands.
                    raise BeyondLogicError(
                        message=(
                            "An antithetical blood group with more than two subtypes has "
                            "reached the point of excluding pairs, and which subtype the "
                            "modifying SNP belongs to is ambiguous. This has no known "
                            "real example - please report the context below."
                        ),
                        context=(
                            f"sample: {bg.sample}, BG: {bg.type}, "
                            f"subtypes: {sorted(d)}, modifying SNP: {modifying_SNP}, "
                            f"zygosity: {bg.variant_pool.get(modifying_SNP)}"
                        ),
                    )
                for pair in bg.alleles[allele_state]:
                    for allele in pair:
                        if allele.number_of_defining_variants == 1:
                            variant = list(allele.defining_variants)[0]
                            if variant.split(":")[1] in antitheticals[bg.type]:
                                to_remove.append(pair)

        if to_remove:
            bg.remove_pairs(
                to_remove, "antithetical_modifying_SNP_is_HOM", allele_state
            )
    return bg


@apply_to_dict_values
def cant_pair_with_ref_cuz_SNPs_must_be_on_other_side(bg: BloodGroup) -> BloodGroup:
    """Filter out allele pairs where reference alleles cannot pair with non-
    reference alleles due to SNP strand requirements.

    This function examines allele pairs from the BloodGroup's NORMAL allele
    state. For each pair that contains a reference allele but not all alleles are
    reference, it splits the pair to analyze the non-reference allele. It then
    determines which SNPs must be on the other strand based on the variant pool's
    zygosity annotations. If a non-reference allele's defining variants are all
    present in the set of SNPs that must be on the other strand (or additional
    criteria are met), the pair is removed.

    Args:
        bg (BloodGroup): A BloodGroup object containing allele states, variant pool,
            and other related data.

    Returns:
        BloodGroup: The updated BloodGroup object after filtering out invalid allele
            pairs.

    Example:
        Given allele pairs such as:
          - JK*01W.03 and JK*01W.04 cannot pair with a reference allele due to SNPs
            like '18:43310313_G_A' and '18:43311054_G_A' requiring specific strand
            orientation.
        Such pairs will be removed from the BloodGroup.

    pair: [Allele(genotype='JK*01W.04',
                  phenotype='.',
                  genotype_alt='.',
                  phenotype_alt='Jk(a+ᵂ)',
                  defining_variants=frozenset({'18:43311054_G_A'}),
                  weight_geno=1000,
                  weight_pheno=3,
                  reference=False,
                  sub_type='JK*01',
                  phases=None,
                  number_of_defining_variants=1),
           Allele(genotype='JK*01', #means HOM for 43319519_ref
                  defining_variants=frozenset({'18:43319519_ref'}),
                  reference=True,

    allele3: Allele(genotype='JK*01W.03',
                    defining_variants=frozenset({'18:43310313_G_A'}),

    bg.variant_pool: {'18:43310313_G_A': 'Heterozygous',
                      '18:43311054_G_A': 'Heterozygous',
                      '18:43311131_G_A': 'Heterozygous',
                      '18:43316538_A_G': 'Heterozygous'}
    SNPS_that_need_to_be_on_other_strand: ['18:43310313_G_A']
    flattened_alleles: {Allele(genotype='JK*01W.03',
                               defining_variants=frozenset({'18:43310313_G_A'}),

                        Allele(genotype='JK*01W.04',
                               defining_variants=frozenset({'18:43311054_G_A'}),

                        Allele(genotype='JK*01N.20',
                               defining_variants=frozenset({'18:43310313_G_A',
                                                            '18:43311054_G_A',
                                                            '18:43311131_G_A',
                                                            '18:43316538_A_G'}),

                        Allele(genotype='JK*01W.11',
                               defining_variants=frozenset({'18:43311054_G_A',
                                                      '18:43310313_G_A'}),
    """
    flattened_alleles = {
        allele
        for allele in flatten_alleles(bg.alleles[AlleleState.NORMAL])
        if not allele.reference
    }

    to_remove = []

    # step1 - for non ref allele find which SNPs have to be on other side
    # (because if they were on same side the allele in question woud be something
    # else)
    for pair in bg.alleles[AlleleState.NORMAL]:
        if pair.contains_reference and not pair.all_reference:
            _, allele = split_pair_by_ref(pair)
            for allele2 in flattened_alleles:
                if allele in allele2:
                    SNPS_that_need_to_be_on_other_strand = [
                        SNP
                        for SNP, zygosity in bg.variant_pool.items()
                        if zygosity == Zygosity.HOM
                    ]

                    het_SNPS = []
                    for SNP in allele2.defining_variants.difference(
                        allele.defining_variants
                    ):
                        if bg.variant_pool[SNP] == Zygosity.HET:
                            het_SNPS.append(SNP)
                    if len(het_SNPS) == 1:
                        SNPS_that_need_to_be_on_other_strand += het_SNPS
                        for SNP in allele.defining_variants:
                            for SNP2 in bg.variant_pool:
                                if SNP == SNP2:
                                    continue
                                if SNP.split("_")[0] in SNP2:
                                    SNPS_that_need_to_be_on_other_strand.append(SNP2)
                        # step2 - if the SNPs that must be on the ref side + and
                        # HOMs define anything then ref not possible
                        for allele3 in flattened_alleles:
                            if all(
                                SNP in SNPS_that_need_to_be_on_other_strand
                                for SNP in allele3.defining_variants
                            ):
                                to_remove.append(pair)
    if to_remove:
        bg.remove_pairs(
            to_remove,
            "cant_pair_with_ref_cuz_SNPs_must_be_on_other_side",
        )

    return bg


@apply_to_dict_values
def ABO_cant_pair_with_ref_cuz_261delG_HET(bg: BloodGroup) -> BloodGroup:
    """Ensure ABO allele pairs obey SNP strand constraints for reference alleles.

    When a reference allele is present (defined by the Lane variant 261del,
    i.e. 9:136132908_T_TC in GRCh37), the paired non-reference allele must be an
    O allele. Some pairings (e.g. ABO*A1.01 with ABO*AEL or ABO*AW variants) are
    not allowed if 136132908_T_TC is present in the defining variants. However,
    pairings with ABO*O alleles are allowed.

    Examples:
        ABO*A1.01/ABO*AEL.02 - Not allowed as 136132908_T_TC is in defining vars.
        ABO*A1.01/ABO*AEL.07 - Not allowed as 136132908_T_TC is in defining vars.
        ABO*A1.01/ABO*AW.25 - Not allowed as 136132908_T_TC is in defining vars.
        ABO*A1.01/ABO*AW.31.01 - Not allowed as 136132908_T_TC is in defining vars.
        ABO*A1.01/ABO*O.01.05 - Allowed.
        ABO*A1.01/ABO*O.01.22 - Allowed.
        ABO*A1.01/ABO*O.01.45 - Allowed.
        ABO*A1.01/ABO*O.01.71 - Allowed.

    Args:
        bg (BloodGroup): A BloodGroup object containing allele pairs, a variant
            pool, and other related genetic data.

    Returns:
        BloodGroup: The updated BloodGroup after filtering out allele pairs where a
            reference allele is improperly paired.
    """
    if bg.type != "ABO":
        return bg
    to_remove = []
    for pair in bg.alleles[AlleleState.NORMAL]:
        if pair.contains_reference and not pair.all_reference:
            ref, allele = split_pair_by_ref(pair)
            # ie, can they exist given other chrom
            if pool_cant_supply_both(ref, allele, bg.variant_pool_numeric):
                to_remove.append(pair)

    if to_remove:
        bg.remove_pairs(
            to_remove,
            "ABO_cant_pair_with_ref_cuz_261delG_HET",
        )

    return bg


def pool_cant_supply_both(ref: Allele, allele: Allele, pool: dict[str, int]) -> bool:
    """Whether one chromosome each is more than the pool has to give.

    Put the reference on one chromosome by spending a copy of each of its defining
    variants, then ask whether the other allele still has a copy of each of its own.
    Where the two alleles share a variant the sample has one copy of, they are asking
    for the same chromosome and the answer is no.

    Only variants the pool holds are spent. A reference allele reaches here straight
    from the database whether or not the sample supports it - NoHomMultiVariantStrategy
    adds it either way - so its variants are routinely absent from a pool built from
    the alleles that were built. Absence is a different claim, and
    cant_name_second_slot_cuz_ref_impossible is the filter that makes it.

    Args:
        ref (Allele): The reference allele of the pair, placed first.
        allele (Allele): The other allele, which must fit on what is left.
        pool (dict[str, int]): Copies per variant. Copied before spending, so the
        caller's dict is left alone.

    Returns:
        bool: True if the two alleles cannot both be present.

    Example:
        HG00436 GYPA. GYPA*02 and GYPA*08 both need 4:144120554_ref, and the sample has
        one copy of it, with 4:144120554_C_A on the other chromosome:

        4:144120554_ref : Heterozygous
        4:144120555_ref : Homozygous
        4:144120567_A_G : Heterozygous
        4:144120567_ref : Heterozygous

        GYPA*02 -> {144120554_ref, 144120555_ref, 144120567_ref}
        GYPA*08 -> {144120554_ref, 144120555_ref, 144120567_A_G}

        Spending GYPA*02 leaves 144120554_ref at 0, so GYPA*08 has nowhere to sit and
        this returns True.
    """
    remaining = dict(pool)
    for variant in ref.defining_variants:
        if variant in remaining:
            remaining[variant] -= 1

    return not all(check_available_variants(0, remaining, operator.gt, allele))


@apply_to_dict_values
def cant_pair_with_ref_cuz_shared_variant_has_too_few_copies(
    bg: BloodGroup,
) -> BloodGroup:
    """Remove a pair whose two alleles want the same chromosome.

    The unphased counterpart of filter_if_all_HET_vars_on_same_side_and_phased, reaching
    the same verdict from copy number instead of phase, so the two arms agree without
    the flag deciding it.

    pair_can_exist is where this arithmetic normally happens, at pair generation, and it
    returns True unchecked for any pair containing the reference. That is not laziness:
    the reference allele arrives from the database rather than from the pool, so its
    variants are often not keys of the pool at all and the subtraction would raise.
    ABO_cant_pair_with_ref_cuz_261delG_HET has done the check properly for one blood
    group for a long time; this is the same check for the rest of them, sharing its
    implementation so the two cannot drift.

    Runs after the ABO filter so ABO pairs keep their own reason string, which is
    documented and user visible. It finds nothing left to do there.

    A '_ref' token counts the same as an alternate and is the larger half of what this
    removes. Heterozygous at a '_ref' token says the reference base is on one chromosome
    and an alternate is on the other, so two alleles both needing it is the same
    contradiction as two alleles both needing one copy of an alternate.

    Args:
        bg (BloodGroup): A BloodGroup object containing allele pairs and a variant pool.

    Returns:
        BloodGroup: The BloodGroup with the impossible pairs removed and recorded.

    Example:
        HG03097 RHCE. RHCE*01 is a reference defined partly by an alternate, because at
        a lane locus the transcript reference differs from the genome reference, and
        RHCE*01.20.01 is RHCE*01 plus one more variant. Both need 1:25420739_G_C:

        1:25390817_G_C : Heterozygous
        1:25390874_ref : Homozygous
        1:25408711_ref : Homozygous
        1:25420739_G_C : Heterozygous
        1:25420739_ref : Heterozygous

        RHCE*01/RHCE*01.20.01 is removed. RHCE*01/RHCE*01.01 and
        RHCE*01.01/RHCE*01.20.01 are untouched, because 1:25420739_ref supplies the
        second chromosome in both.

        HG00619 RHCE is the same thing on a '_ref' token. RHCE*03 and RHCE*01 both need
        1:25408711_ref, which is heterozygous, so the reference base is on one
        chromosome and 1:25408711_G_A is on the other.
    """
    to_remove = []
    for pair in bg.alleles[AlleleState.NORMAL]:
        if pair.contains_reference and not pair.all_reference:
            ref, allele = split_pair_by_ref(pair)
            if pool_cant_supply_both(ref, allele, bg.variant_pool_numeric):
                to_remove.append(pair)

    if to_remove:
        bg.remove_pairs(
            to_remove,
            "cant_pair_with_ref_cuz_shared_variant_has_too_few_copies",
        )

    return bg


def ref_slot_is_impossible(ref: Allele, variant_pool: dict[str, Zygosity]) -> bool:
    """Whether a reference allele needs a variant the sample does not have.

    A reference allele is defined by '_ref' tokens, and one of those is absent from
    the pool when the alternate at that locus is homozygous - there is no reference
    copy left for it to sit on. The two exemptions are the ones no_defining_variant
    already makes, for the same reasons:

    An allele defined only by absence markers ('.') is not making a claim about any
    locus, so there is nothing for the pool to contradict.

    The ABO c.261delG tokens are exempt because the database treats the deletion as
    the reference sequence, which makes ABO*A1.01 a reference allele defined by an
    alternate. Its absence is the ordinary state of a group O sample rather than a
    contradiction.

    Args:
        ref (Allele): The reference allele occupying one slot of a pair.
        variant_pool (dict[str, Zygosity]): The sample's variant pool for this blood
        group. Assumed non-empty - an empty pool is absence of evidence rather than
        evidence of impossibility, and the caller checks that before asking.

    Returns:
        bool: True if the sample cannot carry this reference allele.
    """
    if all(variant.endswith(".") for variant in ref.defining_variants):
        return False

    return any(
        variant not in variant_pool
        for variant in ref.defining_variants
        if variant not in ABO_DELG_VARIANTS
    )


@apply_to_dict_values
def cant_name_second_slot_cuz_ref_impossible(bg: BloodGroup) -> BloodGroup:
    """Name the one slot that is certain when every pair left needs a bad reference.

    The shape: one chromosome matches an allele outright, the other matches nothing
    in the database, and the only pairs the machinery can offer put the reference
    allele in the second slot - which the variant pool says the sample cannot carry.
    Reporting the pair asserts wildtype on a chromosome there is evidence against;
    removing it throws away the slot that *is* known. Neither is what a lab scientist
    writes down. They name the allele they identified and leave the other blank,
    which is UNDETERMINED_SLOT.

    Runs in both arms deliberately. no_defining_variant handles the same
    impossibility but is phased-only (phased.py) and removes the whole pair, so today
    the same sample gets a wrong call unphased and no call phased - a difference
    decided by a flag that is about something else. This runs before it and takes the
    case away from it, so the two arms agree; where other pairs survive this does
    nothing and no_defining_variant behaves as it always has.

    The gate is that *every* remaining pair is this shape. If any pair is callable
    the blood group has an answer and this must not touch it, so a pair with two
    non-reference alleles, a pair of two references, or a reference the pool does not
    contradict all stop it dead. That keeps the change to blood groups whose only
    alternative was a wrong call or no call at all.

    The pair is removed and recorded in filtered_out under this filter's name like
    any other exclusion, so the phenotype engine finds no alleles and reports nothing
    for the blood group. That is deliberate: a phenotype needs both chromosomes, and
    one of them is exactly what the tool has just declined to name.

    Args:
        bg (BloodGroup): The BloodGroup object containing allele pairs and variant pool.

    Returns:
        BloodGroup: The BloodGroup with single_slot_genotypes set and the pairs it came
        from removed, or unchanged if any pair was callable.

    Example:
        GYPA is Lane, so the genome reference is GYPA*02 and 144120567_ref is one of its
        defining variants. A sample called 0|1 at 144120554, 0|1 at 144120555 with the
        caller doubting it, and 1/1 at 144120567:

        FILTER drops 144120555_T_C, and with it GYPA*01, which needs it.

        bg.variant_pool: {'4:144120554_ref': 'Heterozygous',
                          '4:144120555_ref': 'Homozygous',
                          '4:144120567_A_G': 'Homozygous'}

        alleles[NORMAL] -> [Pair(GYPA*02, GYPA*08)]

        4:144120567_ref is not in the pool because the alternate there is homozygous,
        so GYPA*02 has nowhere to sit, while GYPA*08 has all of its defining variants.

        single_slot_genotypes -> ['GYPA*08/Undetermined']
    """
    # .get, not [], because this runs before the co-existing stages create the key. The
    # check is kept anyway so moving the filter later cannot silently start overriding a
    # Knops result.
    if bg.alleles.get(AlleleState.CO) is not None:
        return bg

    pairs = bg.alleles[AlleleState.NORMAL]
    if not pairs or not bg.variant_pool:
        return bg

    named: list[str] = []
    for pair in pairs:
        if not pair.contains_reference or pair.all_reference:
            return bg
        ref, allele = split_pair_by_ref(pair)
        if not ref_slot_is_impossible(ref, bg.variant_pool):
            return bg
        named.append(allele.genotype)

    bg.remove_pairs(
        list(pairs),
        "cant_name_second_slot_cuz_ref_impossible",
        reverts_to_reference=False,
    )
    bg.single_slot_genotypes = [
        f"{genotype}/{UNDETERMINED_SLOT}" for genotype in sorted(set(named))
    ]

    return bg


@apply_to_dict_values
def cant_name_second_slot_cuz_shared_variant_has_too_few_copies(
    bg: BloodGroup,
) -> BloodGroup:
    """Name the slots that are possible when every pair left wants one chromosome twice.

    The sibling of cant_name_second_slot_cuz_ref_impossible, for the other way a pair
    containing the reference can be impossible. There the reference needs a variant the
    sample does not have at all, and only its partner can be named. Here both alleles
    are individually fine and it is the pair that cannot exist, because they share a
    variant the sample has one copy of - so one of them is on that chromosome and the
    other chromosome carries something the database cannot name.

    Which one is on it needs phase to say. Without phase both are candidates, and both
    are reported: naming neither throws away everything that is known, and naming one
    would be a guess. Two candidates contain whatever the phased arm resolves it to, so
    the arms narrow rather than contradict each other.

    Runs immediately before cant_pair_with_ref_cuz_shared_variant_has_too_few_copies and
    on the same predicate, so the two see the same pairs and this takes the case away
    from the other. Where any pair is callable this does nothing and the other filter
    removes the impossible ones as usual.

    The reference is offered as a candidate only where the pool does not contradict it
    on its own. A reference needing a variant that is absent entirely is a different
    claim, and ref_slot_is_impossible is where it is made.

    Args:
        bg (BloodGroup): The BloodGroup object containing allele pairs and variant pool.

    Returns:
        BloodGroup: The BloodGroup with single_slot_genotypes set and the pairs it came
        from removed, or unchanged if any pair was callable.

    Example:
        HG00436 GYPA, long read. The caller emitted no row at 144120555, so GYPA*01 -
        which needs 144120555_T_C - was never built and the only pair left is the one
        that cannot exist:

        4:144120554_ref : c.72G : Heterozygous
        4:144120555_ref : c.71A : Homozygous
        4:144120567_A_G : c.59T>C : Heterozygous
        4:144120567_ref : c.59T : Heterozygous

        Vars_unused:
        4:144120554_C_A : c.72G>T : Heterozygous

        alleles[NORMAL] -> [Pair(GYPA*02, GYPA*08)]

        Both need 4:144120554_ref and the sample has one copy, with 4:144120554_C_A on
        the other chromosome, which completes nothing.

        single_slot_genotypes -> ['GYPA*02/Undetermined', 'GYPA*08/Undetermined']

        The phased arm of the same input reports GYPA*02/Undetermined, which these two
        contain.
    """
    # .get, not [], because this runs before the co-existing stages create the key. The
    # check is kept anyway so moving the filter later cannot silently start overriding a
    # Knops result.
    if bg.alleles.get(AlleleState.CO) is not None:
        return bg

    pairs = bg.alleles[AlleleState.NORMAL]
    if not pairs or not bg.variant_pool:
        return bg

    named: list[str] = []
    for pair in pairs:
        if not pair.contains_reference or pair.all_reference:
            return bg
        ref, allele = split_pair_by_ref(pair)
        if not pool_cant_supply_both(ref, allele, bg.variant_pool_numeric):
            return bg
        named.append(allele.genotype)
        if not ref_slot_is_impossible(ref, bg.variant_pool):
            named.append(ref.genotype)

    bg.remove_pairs(
        list(pairs),
        "cant_name_second_slot_cuz_shared_variant_has_too_few_copies",
        reverts_to_reference=False,
    )
    bg.single_slot_genotypes = [
        f"{genotype}/{UNDETERMINED_SLOT}" for genotype in sorted(set(named))
    ]

    return bg


@apply_to_dict_values
def cant_pair_with_ref_cuz_trumped(bg: BloodGroup) -> BloodGroup:
    """Filter out allele pairs where a reference allele is trumped by a superior allele.

    This function checks allele pairs in the NORMAL state that contain a reference
    allele. If a non-reference allele in the pair is trumped (i.e. has a higher
    weight_geno compared to another allele with one HET variant in the same subtype),
    the pair is removed.


    Args:
        bg (BloodGroup): A BloodGroup object containing allele pairs, variant pool,
            and related genetic data.

    Returns:
        BloodGroup: The updated BloodGroup after filtering out trumped allele pairs.


    Examples:
    ----------
    summary:
        BG Name: FUT3
        Pairs:
            - FUT3*01N.01.02/FUT3*01N.01.12
            - FUT3*01.16/FUT3*01N.01.12
            - FUT3*01/FUT3*01N.01.12
            - FUT3*01/FUT3*01N.01.02
        Filtered Out:
            - FUT3*01/FUT3*01.16 is removed because FUT3*01.16 is trumped by a better allele,
              e.g. FUT3*01N.01.02 or FUT3*01N.01.12.

    BG Name: FUT3
    Pairs:
    Pair(Genotype: FUT3*01N.01.02/FUT3*01N.01.12 Phenotype: ./.)
    Pair(Genotype: FUT3*01.16/FUT3*01N.01.12 Phenotype: ./.)
    Pair(Genotype: FUT3*01/FUT3*01N.01.12 Phenotype: ./.)
    Pair(Genotype: FUT3*01/FUT3*01N.01.02 Phenotype: ./.)

    Filtered Out: defaultdict(<class 'list'>, {'cant_pair_with_ref_cuz_trumped':
    [Pair(Genotype: FUT3*01/FUT3*01.16 Phenotype: ./.)]})

    FUT3*01.16 not with ref due to FUT3*01N.01.02 and FUT3*01N.01.12

    Current truth: FUT3 126

    FUT3*01.16/FUT3*01N.01.12
    FUT3*01/FUT3*01N.01.02
    FUT3*01/FUT3*01N.01.12
    FUT3*01N.01.02/FUT3*01N.01.12
    [Allele(genotype='FUT3*01.15',
            defining_variants=frozenset({'19:5844184_C_T',
                        '19:5844367_C_T'}), HOM
            weight_geno=1000,
            reference=False,
    Allele(genotype='FUT3*01.16',
            defining_variants=frozenset({'19:5844043_C_T',
                                        '19:5844184_C_T',
                                        '19:5844367_C_T'}), HOM
            weight_geno=1000,
            reference=False,
    Allele(genotype='FUT3*01N.01.02',
            defining_variants=frozenset({'19:5844184_C_T',
                                        '19:5844367_C_T', HOM
                                        '19:5844838_C_T'}), HOM
            weight_geno=1,
            reference=False,
    Allele(genotype='FUT3*01N.01.12',
            defining_variants=frozenset({'19:5843883_C_G',
                                        '19:5844367_C_T', HOM
                                        '19:5844838_C_T'}), HOM
            weight_geno=1,
            reference=False,

    res[goi].variant_pool: {'19:5843883_C_G': 'Heterozygous',
                            '19:5844043_C_T': 'Heterozygous',
                            '19:5844184_C_T': 'Heterozygous',
                            '19:5844367_C_T': 'Homozygous',
                            '19:5844838_C_T': 'Homozygous'}
    """

    to_remove = []
    flattened_alleles = flatten_alleles(bg.alleles[AlleleState.NORMAL])

    if not any(allele.reference for allele in flattened_alleles):
        return bg
    alleles_without_ref = [
        allele for allele in flattened_alleles if not allele.reference
    ]

    alleles_with_1_HET_var = [
        allele
        for allele in alleles_without_ref
        if [bg.variant_pool[variant] for variant in allele.defining_variants].count(
            Zygosity.HET
        )
        == 1
    ]
    if alleles_with_1_HET_var:
        for pair in bg.alleles[AlleleState.NORMAL]:
            if pair.contains_reference and not pair.all_reference:
                ref, allele = split_pair_by_ref(pair)
                for allele2 in alleles_with_1_HET_var:
                    if allele.sub_type != allele2.sub_type:
                        continue
                    if allele.weight_geno > allele2.weight_geno:
                        to_remove.append(pair)
    if to_remove:
        bg.remove_pairs(to_remove, "cant_pair_with_ref_cuz_trumped")

    return bg


@apply_to_dict_values
def cant_not_include_null(bg: BloodGroup) -> BloodGroup:
    """Remove allele pairs without null alleles when a null with ≤1 HET variant exists.

    Filters out non-null allele pairs when there exists at least one null allele
    that has one or fewer heterozygous variants. This ensures that likely null
    alleles (those with minimal heterozygous support) are prioritized in the
    final allele configuration.

    Args:
        bg: BloodGroup object containing allele pairs and variant information.

    Returns:
        The modified BloodGroup object with non-null pairs removed if applicable,
        or the original object unchanged if no qualifying null alleles exist.

    Note:
        Only processes alleles in the NORMAL state. Pairs are removed via the
        remove_pairs method with reason "cant_not_include_null" for tracking.
    """
    flattened_alleles = flatten_alleles(bg.alleles[AlleleState.NORMAL])
    null_alleles = [allele for allele in flattened_alleles if allele.null]
    null_alleles_with_one_HET_variant = [
        allele for allele in null_alleles if one_HET_variant(allele, bg.variant_pool)
    ]
    if not null_alleles_with_one_HET_variant:
        return bg
    to_remove = []
    for pair in bg.alleles[AlleleState.NORMAL]:
        if pair.allele1.null or pair.allele2.null:
            continue
        to_remove.append(pair)
    if to_remove:
        bg.remove_pairs(to_remove, "cant_not_include_null")

    return bg


@apply_to_dict_values
def cant_have_2_non_ref_alleles_cuz_only_1_gene_copy(bg: BloodGroup) -> BloodGroup:
    """Two different non-reference alleles need two gene copies to sit on.

    Where a caller reported one copy of the gene on two chromosomes, one chromosome
    carries the gene and the other carries no gene at all. So there is exactly one place
    for a real allele, and a pair naming two different ones is claiming both are on it.

    Before this filter it was noticed at the reporting layer and warned about, then the
    pair was written out as called anyway - the one impossible result the pipeline
    knowingly reported. A warning is not an audit trail: it does not say which pair, it
    is not in filtered_out, and it does not reach the debug trace or the PDF. Now the
    pair is excluded by name like every other impossible pair.

    Two things it deliberately leaves alone:

    - **Co-existing alleles.** Two alleles on one chromosome is exactly what co-existing
      *means*, so on the Knops path the shape is legitimate rather than impossible, and
      it is what gets reported when it is present. Skipped entirely rather than filtered
      and then ignored, so it cannot empty NORMAL and trigger the reverted-to-reference
      warning for a blood group whose answer was never coming from NORMAL.
    - **The same allele twice.** The pairing machinery writes a duplicate where one copy
      carries an allele and there is nothing to pair it with, and get_genotypes collapses
      that to 'allele/absent_slot'. One allele, one copy, no contradiction - which is why
      this counts *distinct* genotypes rather than alleles.

    Nothing here fires unless a caller reported gene copy number, so on an input without
    that channel locus_copies is None and this is a no-op.

    Args:
        bg (BloodGroup): The BloodGroup object containing allele pairs.

    Returns:
        BloodGroup: The updated BloodGroup, with impossible pairs moved to filtered_out
        under 'cant_have_2_non_ref_alleles_cuz_only_1_gene_copy'.
    """
    if bg.locus_copies != 1 or bg.chrom_copies != 2:
        return bg
    # .get, not [], on both: bg.alleles is a plain dict and the states it holds depend on
    # how far down the pipeline the blood group is.
    if bg.alleles.get(AlleleState.CO) is not None:
        return bg
    if not bg.alleles.get(AlleleState.NORMAL):
        return bg

    to_remove = []
    for pair in bg.alleles[AlleleState.NORMAL]:
        non_ref = {
            allele.genotype for allele in pair.alleles if not allele.reference
        }
        if len(non_ref) > 1:
            to_remove.append(pair)
    if to_remove:
        bg.remove_pairs(
            to_remove, "cant_have_2_non_ref_alleles_cuz_only_1_gene_copy"
        )

    return bg


@apply_to_dict_values
def ensure_HET_SNP_used(bg: BloodGroup) -> BloodGroup:
    """
    Ensures that heterozygous variants are utilized in allele pairs if they can form
    existing alleles.

    This function iterates over heterozygous variants in the variant pool of a
    BloodGroup object.
    For each heterozygous variant, it checks whether adding this variant to the defining
    variants of alleles in each pair results in alleles that already exist outside the
    pair. If multiple such matches are found, the pair is considered invalid and is
    removed from the allele pairs.

    Args:
        bg (BloodGroup): The BloodGroup object containing allele pairs and variant pool
        information.

    Returns:
        BloodGroup: The updated BloodGroup object with inconsistent allele pairs
        removed.

    Example:
        Suppose you have a BloodGroup object with allele pairs and a variant pool that
        includes heterozygous variants. This function will remove pairs where the
        heterozygous variants are not properly used in the allele definitions, but
        could form existing alleles when added.


        #Results:
        Genotypes count: 2
        Genotypes:
        FUT3*01N.03.02/FUT3*01N.03.02 #not possible as 19:5843969_G_A isn't
        accounted for
        FUT3*01N.03.02/FUT3*01N.03.09
        Phenotypes (numeric):
        Phenotypes (alphanumeric): Le(a-b-)

        #Data:
        Vars:
        19:5844537_ref : Homozygous
        19:5844649_ref : Homozygous
        19:5843969_G_A : Heterozygous
        Vars_phase:

        Vars_phase_set:

        Raw:
        Allele
        genotype: FUT3*01.04
        defining_variants:
                19:5844537_ref
        weight_geno: 1000
        phenotype: . or Active
        reference: False

        Allele
        genotype: FUT3*01N.03.01
        defining_variants:
                19:5844649_ref
        weight_geno: 1
        phenotype: . or Le(a-),Le(b-)
        reference: False

        Allele
        genotype: FUT3*01N.03.02
        defining_variants:
                19:5844649_ref
                19:5844537_ref
        weight_geno: 1
        phenotype: . or Le(a-),Le(b-)
        reference: False

        Allele
        genotype: FUT3*01N.03.09
        defining_variants:
                19:5844649_ref
                19:5843969_G_A
                19:5844537_ref
        weight_geno: 1
        phenotype: . or Le(a-),Le(b-)
        reference: False
    """
    for allele_state in [AlleleState.NORMAL, AlleleState.CO]:
        if not proceed(bg, allele_state):
            continue
        to_remove = []
        for variant, zygo in bg.variant_pool.items():
            if zygo == Zygosity.HET:
                for pair in bg.alleles[allele_state]:
                    if all_hom(bg.variant_pool, pair.allele1) and all_hom(
                        bg.variant_pool, pair.allele2
                    ):
                        hits = check_var(bg, pair, allele_state, variant)
                        if hits:
                            to_remove.append(pair)
        if to_remove:
            bg.remove_pairs(to_remove, "ensure_HET_SNP_used", allele_state)

    return bg


@apply_to_dict_values
def filter_HET_pairs_by_weight(bg: BloodGroup) -> BloodGroup:
    """This filter forced us to make decisions where there was not always
    a clearly correct answer. I have left the example long to show the evolution of
    thought and why we landed where we did.

    Args:
        bg (BloodGroup): The BloodGroup object containing allele pairs and variant pool
        information.

    Returns:
        BloodGroup: The updated BloodGroup object with inconsistent allele pairs
        removed.

    Example
    ----------
     FUT1/2
    All bar ref here is HET
    So while it's possible that any allele can exist in a pair,
    no pair can exist without most weighted allele, if it just has 1 defining SNP
    (or should it be just has 1 HET defining SNP, cuz if the rest HOM, same, right?)
    # update - only if alleles are same subtype! and if the 1 SNP of the most weighted
    # allele is 'in' one of the pair (meaning that most weighted allele doesn't exist?)
    . Ref has to be included (as all SNPS could be together)
    but FUT2*01/FUT2*01.03.01 isn't possible
    Further, FUT2*01/FUT2*01N.16 isn't possible as FUT2*01.03.01 trumps FUT2*01N.16
    [Allele(genotype='FUT2*01',
        genotype_alt='.',
        defining_variants=frozenset({'19:49206250_ref'}), hom
        weight_geno=LOW_WEIGHT,
        weight_pheno=2,
        reference=True,
        Sub_type='FUT2*01'),
    Allele(genotype='FUT2*01.03.01',
        genotype_alt='.',
        defining_variants=frozenset({'19:49206286_A_G'}),
        weight_geno=LOW_WEIGHT,
        weight_pheno=1,
        reference=False,
        Sub_type='FUT2*01'),
    Allele(genotype='FUT2*01N.02',
        genotype_alt='.',
        defining_variants=frozenset({'19:49206674_G_A'}),
        weight_geno=1,
        weight_pheno=5,
        reference=False,
        Sub_type='FUT2*01'),
    Allele(genotype='FUT2*01N.16',
        genotype_alt='.',
        defining_variants=frozenset({'19:49206985_G_A'}),
        weight_geno=8,
        weight_pheno=5,
        reference=False,
        Sub_type='FUT2*01')]

    res["FUT2"].variant_pool: {'19:49206250_ref': 'Homozygous',
                               '19:49206286_A_G': 'Heterozygous',
                               '19:49206674_G_A': 'Heterozygous',
                               '19:49206985_G_A': 'Heterozygous'}


    issue

    [Allele(genotype='JK*01N.19',
                        genotype_alt='.',
                        defining_variants=frozenset({'18:43319274_G_A'}),
                        weight_geno=7,
                        reference=False,
                        sub_type='JK*01',

    Allele(genotype='JK*02N.17',
                        genotype_alt='.',
                        defining_variants=frozenset({'18:43319274_G_A',
                                                    '18:43319519_G_A'}),
                        weight_geno=13,
                        reference=False,
                        sub_type='JK*02',

    Allele(genotype='JK*01W.06',
                        genotype_alt='.',
                        defining_variants=frozenset({'18:43310415_G_A',
                                                    '18:43316538_A_G'}), hom
                        weight_geno=1000,
                        reference=False,
                        sub_type='JK*01',

    JK*02	43319519_G_A

    pool: {'18:43310415_G_A': 'Heterozygous',
                '18:43316538_A_G': 'Homozygous',
                '18:43319274_G_A': 'Heterozygous',
                '18:43319519_G_A': 'Heterozygous'}

    current code removed JK*01W.06/JK*02N.17. At first glance I thought this was the
    correct behaviour
    because JK*01N.19 trumps JK*01W.06 and only has 1 SNP. but I noticed
    that the 1 defining SNP for JK*01N.19 (43319274_G_A) is part of the definition of
    JK*02N.17, which is further complicated by the fact its a different subtype.
    Based on my understanding I think that JK*01W.06/JK*02N.17 is right, because
    the trumpiness of JK*01N.19 is lost to JK*02N.17.

    Is there a clear rule to cover this problem? One the one hand, nulls are always
    the most weighted allele (they trump other options) but at the same time alleles
    whose defining variants are all in the list of defining variants of another alleles
    are no longer considered. Which of these rules is applied first? In this example
    the 1 defining SNP for JK*01N.19 (43319274_G_A) is part of the definition of
    JK*02N.17, so JK*01N.19 'doesn't exist' anymore. That makes sense to me, but if
    JK*01N.19 caused truncation or frameshift then I would say that this should
    be the most important consideration... unless the addition of the other variant to
    make JK*02N.17 rescues from being null? As I think this through I think I'm flipping
    and would think that any SNP that makes null would have to be to first consideration

    assumig we want to handle 'in' and subtype:
        each pair will have just 1 of each subtype (other filters will remove)
        so if pair has to include most weighted allele of subtype (sigh)
    """

    flattened_alleles = flatten_alleles(bg.alleles[AlleleState.NORMAL])
    if not flattened_alleles:
        return bg
    if all(allele.weight_geno == LOW_WEIGHT for allele in flattened_alleles):
        return bg
    sub_types = {allele.sub_type for allele in flattened_alleles}

    keepers = []

    for sub_type in sub_types:
        max_weight = min(
            [
                allele.weight_geno
                for allele in flattened_alleles
                if allele.sub_type == sub_type
            ]
        )  # lowest weight, base 1

        weights_with_1_SNP = [
            allele.weight_geno
            for allele in flattened_alleles
            if one_HET_variant(allele, bg.variant_pool) and allele.sub_type == sub_type
        ]
        max_weight_with_1_SNP = (
            min(weights_with_1_SNP) if weights_with_1_SNP else LOW_WEIGHT
        )

        if max_weight_with_1_SNP == max_weight:
            alleles_with_max_weight_and_1_SNP = [
                allele
                for allele in flattened_alleles
                if one_HET_variant(allele, bg.variant_pool)
                and allele.sub_type == sub_type
                and allele.weight_geno == max_weight
            ]
            for pair in bg.alleles[AlleleState.NORMAL]:
                if all(
                    allele in pair.allele1 or allele in pair.allele2
                    for allele in alleles_with_max_weight_and_1_SNP
                ) or any(
                    allele.weight_geno == max_weight
                    for allele in pair
                    if allele.sub_type == sub_type
                ):
                    keepers.append(pair)
        else:  # if the allele with max weight does't have just 1 het SNP, need to keep all
            for pair in bg.alleles[AlleleState.NORMAL]:
                keepers.append(pair)
    to_remove = [pair for pair in bg.alleles[AlleleState.NORMAL] if pair not in keepers]

    if to_remove:
        bg.remove_pairs(to_remove, "filter_HET_pairs_by_weight")

    return bg


@apply_to_dict_values
def filter_pairs_by_context(bg: BloodGroup) -> BloodGroup:
    """Filter allele pairs by context.

    This function removes allele pairs from the NORMAL state of the BloodGroup if
    the context indicates that the pair cannot exist. It checks whether the
    combination of remaining variant counts and an allele's defining variants
    match those in other alleles. This ensures that pairs like
    'A4GALT*01/A4GALT*02' are filtered out when a more comprehensive allele is
    implied by the available variants. This logic overlaps with
    cant_pair_with_ref_cuz_SNPs_must_be_on_other_side, but is still needed

    Args:
        bg (BloodGroup): A BloodGroup object containing allele pairs and a numeric
            variant pool.

    Returns:
        BloodGroup: The updated BloodGroup with contextually invalid pairs removed.

    Example:
        Given allele definitions:
            - A4GALT*01.02: defining_variants = {'22:43089849_T_C'}
            - A4GALT*02:   defining_variants = {'22:43113793_C_A'}
            - A4GALT*02.02: defining_variants = {'22:43113793_C_A',
                                                 '22:43089849_T_C'}
            - A4GALT*01:   defining_variants = {'22:43113793_ref'}
        And a variant pool:
            {'22:43089849_T_C': 'Heterozygous',
             '22:43113793_C_A': 'Heterozygous'}
        The valid pairs are:
            'A4GALT*01.02/A4GALT*02' and
            'A4GALT*01/A4GALT*02.02',
        while 'A4GALT*01/A4GALT*02' is not possible.
    """
    # need to catch; if remaing var and allele var combine to define another allele
    # and if remaining var defines or combines with other allele to define another
    # allele then pair is not possible, but check that the pair that rules it out
    # does already exist

    to_remove = []
    for pair in bg.alleles[AlleleState.NORMAL]:
        def_vars = {
            a.defining_variants
            for pair2 in bg.alleles[AlleleState.NORMAL]
            for a in pair2
            if a not in pair
        }
        for allele in pair:
            variant_pool_copy = bg.variant_pool_numeric.copy()
            if allele.reference:
                continue

            for variant in allele.defining_variants:
                variant_pool_copy[variant] -= 1
            left_over_vars = [k for k, v in variant_pool_copy.items() if v > 0]
            if any([len(left_over_vars) == 0, len(def_vars) < 2]):
                continue
            remaining = [
                tuple(sorted(set(left_over_vars))),
                tuple(sorted(set(left_over_vars + list(allele.defining_variants)))),
            ]
            if all(variants in def_vars for variants in remaining):
                to_remove.append(pair)
                break
    if to_remove:
        bg.remove_pairs(to_remove, "filter_pairs_by_context")

    return bg


@apply_to_dict_values
def impossible_alleles(bg: BloodGroup) -> BloodGroup:
    """
    Rule out impossible alleles based on Homozygous variants.

    Args:
        bg (BloodGroup): A BloodGroup object containing allele states and phasing
            information.

    Returns:
        bg (BloodGroup): A BloodGroup object containing allele states and phasing
            information.

    Example:

    Allele(genotype='JK*02W.03',
            defining_variants={'18:43310415_G_A',
                                18:43316538_A_G',
                               '18:43319519_G_A'},
    Allele(genotype='JK*02W.04',
           defining_variants={'18:43310415_G_A',
                              '18:43319519_G_A'},

    '18:43310415_G_A': 'Heterozygous',
    '18:43316538_A_G': 'Homozygous',
    '18:43319519_G_A': 'Heterozygous',

    JK*02W.04 is impossible because '18:43316538_A_G' is Homozygous
    """
    for allele_state in [AlleleState.NORMAL, AlleleState.CO]:
        if not proceed(bg, allele_state):
            continue
        alleles = list(flatten_alleles(bg.alleles[allele_state]))
        homs = {
            variant for variant, zygo in bg.variant_pool.items() if zygo == Zygosity.HOM
        }
        impossible_alleles = []
        for allele in alleles:
            for allele2 in alleles:
                if allele.genotype != allele2.genotype and allele2 in allele:
                    dif = allele.defining_variants.difference(allele2.defining_variants)
                    if all(variant in homs for variant in dif):
                        impossible_alleles.append(allele2)

        to_remove = []
        for pair in bg.alleles[allele_state]:
            if pair.allele1 in impossible_alleles or pair.allele2 in impossible_alleles:
                to_remove.append(pair)
        if to_remove:
            bg.remove_pairs(to_remove, "filter_impossible_alleles", allele_state)
    return bg
