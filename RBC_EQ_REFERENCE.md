# RBCeq2 implementation reference

The body below was moved verbatim from an AGENTS.md snapshot saved on 2026-09-07.
Its implementation/status text predates the correction in commit `fe6a9a7331e1365dea8e4c8f91fbad0cbe3e1fa0`.
Read only the sections relevant to the task; current operating rules are authoritative in AGENTS.md.
Implementation/status claims, measurements, and source line numbers below are historical:
verify current code before relying on them. This date note does not retire the applicable constraints.

Current correction: the old "no reference equivalent" statement in the `pair_can_exist` discussion
predates the filters `cant_pair_deletion_with_ref_cuz_HEM_SNP_defines_an_allele` and
`cant_pair_with_ref_cuz_a_deletion_names_the_missing_copy`. Inspect those functions in current source.

## Topic index

| Read for | Line in this file |
|---|---|
| [Repo map](#repo-map) | 35 |
| [Architecture](#architecture): pipeline order and stage shape | 63 |
| [Domain glossary](#domain-glossary) | 90 |
| [Code conventions](#code-conventions) | 110 |
| [Known landmines](#known-landmines): read those related to touched code | 123 |
| Haploid GT, locus/chromosome/token copies, slot markers, dosage, no-calls | 127 |
| Reference pairs, token spending, named exclusion/rescue filters | 187 |
| Missing tokens, locus evidence, `NO_COPIES` versus `NO_DATA` | 243 |
| `FILTER` meanings, exposure, and classification table | 264 |
| Assertions, antithetical subtypes, DB whitespace, and test doubles | 286 |
| Multi-allelic encoding, duplicate rows, and plural `FILTER` lookup | 321 |
| Structural token coordinates, event/definition collisions, RH and SV scope | 361 |
| Output delimiters, alternatives, and cross-file ordering | 384 |
| [Design constraints inherited from history](#design-constraints-inherited-from-history) | 418 |
| [When changing filter logic](#when-changing-filter-logic) | 434 |

<!-- Original AGENTS.md section content begins below; retained byte for byte. -->

## Repo map

```
src/rbceq2/
  main.py                    entry point; composes the whole pipeline as a list of partials
  core_logic/
    data_procesing.py        NOTE THE TYPO (one 's'). Variant pool, zygosity, pair generation
    alleles.py               Allele, Line, BloodGroup, Pair dataclasses
    co_existing.py           Knops-specific co-existing allele logic
    large_variants.py        SV reading (SvReader) and fuzzy matching (SvMatcher)
    constants.py             AlleleState, LANE, antithetical maps, BG name mappings
    utils.py                 Zygosity, decorators, BeyondLogicError
  filters/
    geno.py                  unphased genotype filters
    phased.py                phased equivalents
    knops.py                 co-existing filters
  phenotype/
    antigens.py              antigen objects, modifiers, antithetical reconciliation
    choose_pheno.py          genotype pair -> phenotype
  IO/
    vcf.py                   VCF reading, phase sets, lane variants (native gVCF unsupported)
    record_data.py           TSV output
    PDF_reports.py           optional per-sample PDF
  db/db.py                   DB loading, reference alleles, lane variants
  resources/db.tsv           the allele database
tests/                       one test_*.py per module
```

## Architecture

`main.py` builds an ordered list of callables and `reduce`s them over a dict of
`{BG name: BloodGroup}`. Order matters and is load-bearing — one filter carries the comment
`has to be after normal filters!!!!!!!`. Broad shape:

1. **Read** — VCF parsed; ploidy inferred per sample (`VCF._infer_haploid_chroms`, before
   `remove_home_ref` because it decides whether a haploid `0` is a hom ref call); SVs extracted
   (`SvReader`) and fuzzy-matched to DB tokens (`SvMatcher`)
2. **Raw alleles** — every DB allele whose defining variants are all present
3. **Variant pool** — `make_variant_pool` assigns zygosity per variant token, and sets
   `bg.chrom_copies` from `Db.single_copy_types` × `VCF.haploid_chroms`
4. **Large-indel adjustment** — `modify_variant_pool_if_large_indel` rewrites HOM→HEM inside deletions
5. **Phasing** — phase pool built from raw GT strings; unphased alleles removed if `--phased`
6. **Pairs** — all viable pairs generated from the pool
7. **Filters** — ~23 named filters, plus Knops co-existing logic
8. **Phenotype** — antigen merge, modifier weighting, antithetical reconciliation
9. **Output** — 3 TSVs (genotype, phenotype numeric, phenotype alphanumeric), log, optional PDFs

Pairs stay two-slot all the way through, including for a single-copy region, where both slots
hold the same allele. The haploid shape is a *reporting* decision and appears only at
`get_genotypes` and `add_refs`, which write `HAPLOID_SECOND_SLOT` (`-`) into the second slot —
`XK*N.03/-`. Nothing in the filters or the phenotype engine sees it.

Most stage functions are decorated with `@apply_to_dict_values`, which maps the function over
every BloodGroup. Write new stages the same way.

## Domain glossary

| term | meaning |
|---|---|
| BG | blood group |
| ISBT | International Society for Blood Transfusion — defines the allele nomenclature |
| SNV / SNP | function names use `SNP`, prose uses `SNV`. Same thing |
| `_ref` | a defining variant meaning "no change from reference at this position" |
| antithetical antigens | two antigens that are opposites (K vs k). Homozygosity for one excludes the other |
| trumps | within a subtype, a higher-weighted allele displaces a lower one on the same chromosome |
| "allele A is *in* allele B" | A's defining variants are a subset of B's |
| co-existing | Knops-only: two alleles combining on one chromosome, written with `+` |
| weight | genotype rank; **1 is highest**, 1000 is the default/lowest. Nulls outrank non-nulls |
| HOM / HET / HEM | homozygous / heterozygous / hemizygous |
| lane variant | a locus where the transcript reference differs from the genome reference |
| `chrom_copies` | copies of a *region* the sample was born with → **allele slots in the result**. 2, or 1 for non-PAR X/Y in a male. Never changed by a deletion |
| `locus_copies` | copies of the sequence still present at a locus, ie after a deletion. 0–2. This is what `HEM` (1) and `NO_COPIES` (0) encode |
| `token_copies` | how many of those carry *this* token. Zero is encoded as absence from the pool |
| PAR | pseudoautosomal — X/Y regions present twice in everyone. XG and CD99 are inside PAR1 and stay diploid; XK, GATA1 and ATP11C are outside it |

## Code conventions

- **Abbreviations in function names are intentional.** `cuz` = because. Filter names are long
  and descriptive on purpose: `cant_pair_with_ref_cuz_SNPs_must_be_on_other_side`. Match the
  style; do not "clean up" existing names — they appear verbatim in user-facing debug logs and
  in the published documentation.
- Variant token format: `chrom:pos_REF_ALT`, e.g. `18:43319519_C_T`; `_ref` for no change;
  large variants as `pos_del_59kb` (rounded to the nearest kb deliberately — exact sizes are
  too brittle given breakpoint wobble).
- `chr` prefixes are stripped on read; chromosomes are `1`–`22`, `X`, `Y`.
- Logging via `loguru`. Warnings are surfaced to users deliberately — they usually indicate an
  input problem or something novel. Do not downgrade warnings to debug to quiet output.

## Known landmines

Things that look fine and are not. Verify against these before touching related code.

- **A haploid GT has three readings, and two of them are live.** The same `1` means different
  things in different places and the results have **different shapes**, so `get_ref`
  (`data_procesing.py:1363`) decides by which of two counts is 1 rather than by the GT:

  | which count is 1 | what it means | allele slots | second slot |
  |---|---|---|---|
  | `chrom_copies` | outside PAR on X/Y; the sample has one chromosome | **1** | `-` (`HAPLOID_SECOND_SLOT`) |
  | `locus_copies` | the caller encoded gene copy number as GT ploidy | **2** | the DB's subtype for a missing copy, else `Novel_gene_deletion` (`NOVEL_DELETION_SLOT`) |
  | neither | a haploid GT where the sample has two of everything | — | refuses by name, `get_ref/haploid_GT_where_neither_count_is_one` |

  Rows B1 and D2 of `ploidy_state_table.md` are the *same GT* with different correct answers,
  which is the whole reason the two counts are separate. Reading the second as the first
  reports `RHD*01/-` and loses a chromosome the sample has. `locus_copies_for_bg`
  (`:430`) demands agreement across every database locus the VCF reported for that gene, so
  one haploid GT among diploid neighbours stays a caller quirk and is refused, not read as
  copy number.
- **`locus_copies` cannot say zero.** It returns 1 or `None` and nothing else, so a gene with no
  copies at all has no way to be signalled through it. `Zygosity.NO_COPIES` covers that case
  only when a deletion record exists, because `modify_variant_pool_if_large_indel` builds its
  `big_dels` from tokens already in the pool (`data_procesing.py:864`). Input that reports copy
  number without calling structural variants — an array — therefore has no CN 0 path at all,
  which is row D1 and still open. The naming half is solved (`Db.get_gene_absent_subtypes`,
  per gene, a subtype rather than an allele because a copy number carries no breakpoints); the
  detection half is not. Do not assume a CN 0 signal exists to read.
- **The three second-slot markers are three different claims.** `-` says there is no second
  chromosome, `Novel_gene_deletion` says there is one and it carries no copy of the gene,
  `Undetermined` says it
  carries a copy whose allele the database cannot name. Collapsing any pair of them undoes the
  distinction `chrom_copies`/`locus_copies` exists to draw, and pairing with the reference
  instead — the option that looks tidiest — asserts wildtype on a chromosome there is positive
  evidence against. All three are user visible and documented; see `constants.py:584-615`.
- **Ploidy is per genotype, not per file, and above two copies the question is dosage.**
  Nothing in a VCF header declares ploidy; it is the number of allele indices in each `GT`,
  so it varies legitimately between samples and between records. A genotype naming four
  copies is ordinary input, and a gene conversion caller reporting a paralogue pair writes
  one. It is read where its **dosage** — the count of alternate copies — is 0 or equals
  the ploidy, because only then does it project onto the chromosomes the sample has
  whichever copies go where; anything between refuses by name,
  `get_ref/dosage_between_the_bounds` (`zygosity_of_non_diploid_GT`,
  `data_procesing.py:1565`). Three traps in one branch:
  - **Count, never pattern match.** Ascending order in an unphased genotype is a convention,
    not a rule, so `1/0/1/1` and `0/1/1/1` are the same statement and a shape test misses
    one.
  - **The missing genotype scales too** — `./././.` is what `./.` is at four copies, and
    it read as an error until the no-call test was moved ahead of the ploidy gate.
  - **`HOM_REF_GTS` covers the high-ploidy spellings by prefix, not by decision.**
    `'0/0/0/0'.startswith('0/0')` is what drops those rows in `remove_home_ref`
    (`vcf.py:587`); "tidying" that into exact matching reopens the gap silently. The
    membership test in `get_variants` (`:848`) is exact and does *not* cover them, which is
    harmless only because dosage 0 now encodes as absence rather than as a homozygous ALT.
  `Number=G` fields (`PL`, `GL`, `GP`) grow combinatorially with ploidy, which is irrelevant
  here: **only `GT` and `PS` are ever read off a row.**
- **`./.` is overloaded two ways**: a genuine no-call and, from copy-number-aware callers, zero
  copies. These mean opposite things. The synthesised lane row is *not* one of them any more —
  it carries `SYNTHESISED_HOM_REF_GT` (`constants.py:23`), which is the only legitimate
  assertion of wildtype in the codebase.
- **`Zygosity` is a projection of two numbers, not one label.** `HEM` is one copy of the locus
  carrying one copy of the token; `NO_COPIES` is zero of both. Neither says how many allele
  slots the *result* has — that is `BloodGroup.chrom_copies`, a third number, and it is never
  changed by a deletion. Reading pair shape off a per-token zygosity is the bug issue #40 was.
- **`pair_can_exist` never checks a pair containing the reference, and must not be made to.**
  It subtracts one allele's defining variants from the pool and requires a copy left for the
  other's, but returns True unchecked whenever either allele is the reference
  (`data_procesing.py:2502`, short circuit at `:2520-2522`). Deleting those three lines fails
  *every* sample with `KeyError`: the reference arrives from the database rather than from the
  pool — `NoHomMultiVariantStrategy` (`:2361`) adds it whether the sample supports it or not —
  so its tokens are routinely not keys of a pool built from the alleles that were *built*.

  **That `KeyError` is a property of the naive edit, not of the problem, and it is not the
  reason to leave this alone.** Skipping absent tokens on both sides — spending and checking,
  which is what `pool_cant_supply_both` already does — removes the short circuit with no
  `KeyError` anywhere, and rejects 256 reference pairs on 17 samples. **The reason is that a
  pair rejected at generation cannot be rescued.** Where removing every impossible pair would
  leave a blood group with nothing, `cant_name_second_slot_cuz_shared_variant_has_too_few_copies`
  names both candidates instead; `pair_can_exist` sees one candidate pair in isolation, before
  anything knows how many will survive, and can only decline to build it. Measured: moving the
  arithmetic there leaves the truth set and the 967-sample set byte identical — neither
  exercises the rescue — and on the long read set turns `HG01527` RHCE from
  `RHCE*01.20.04.02/Undetermined`, and `HG03761` from two named candidates, into
  `Undetermined/Undetermined`. `HG01527` is the sample item 5 was closed on. Even where
  answers hold, ~1,500 exclusions on the 967-sample set stop being named: 903 from
  `cant_pair_with_ref_cuz_shared_variant_has_too_few_copies` and 606 from
  `filter_pairs_on_antithetical_zygosity`, silently, for identical output. **Covering each gap
  with a named filter is the design, not a workaround** — it is what buys the exclusion trail.

  Two filters cover the gap instead, both spending only the tokens the pool holds:
  `cant_pair_with_ref_cuz_shared_variant_has_too_few_copies` (`geno.py:473`) removes the pair,
  and `cant_name_second_slot_cuz_shared_variant_has_too_few_copies` (`:661`) runs one line
  earlier and names both candidates where removing would leave the blood group with nothing.
  `pool_cant_supply_both` (`:426`) is the arithmetic, lifted out of
  `ABO_cant_pair_with_ref_cuz_261delG_HET` (`:380`), which had been doing this correctly for
  one blood group all along. Five things that are easy to get wrong here:
  - **A `_ref` token constrains exactly as hard as an alternate, and is the larger half.**
    Heterozygous at a `_ref` token says the reference base is on one chromosome and an
    alternate on the other, so two alleles both needing it is the same contradiction as two
    alleles both needing one copy of an alternate. Of the 10,854 genotype strings this removed
    across nine datasets, 6,318 were `_ref` only against 2,861 on an alternate alone.
  - **It is not confined to the four references defined partly by an alternate.** `KN*01`,
    `RHD*01`, `ABO*A1.01` and `RHCE*01` are the only four, on both builds — but GYPA's 2,821
    come from `GYPA*02`, an ordinary all-`_ref` reference.
  - **A token *absent* from the pool is a different claim**, made by `ref_slot_is_impossible`
    (`:539`) and acted on by `no_defining_variant`. `pool_cant_supply_both` skips what the pool
    does not hold rather than reading absence as zero, which is also why it cannot raise. See
    the next landmine for why absence on its own settles nothing.
  - **Reason strings moved for exclusions that were already being made.** Running earlier than
    `filter_impossible_alleles`, `cant_pair_with_ref_cuz_SNPs_must_be_on_other_side` and
    `filter_if_all_HET_vars_on_same_side_and_phased`, the new filter is now credited with
    cases those used to name — 59 in `dragen_per_sample`, 184 in `dragen_joint_3209`. Same
    answers, more specific reason, but the strings are user visible and documented.
  - **All of this is about *token* copies, and says nothing about *gene* copies.** The pool is
    the only thing any of it sees, so none of it notices that a reference has been paired onto
    a chromosome carrying no copy of the gene at all.
    `cant_have_2_non_ref_alleles_cuz_only_1_gene_copy` (`geno.py:901`) is the gene copy check
    and covers only the two-non-reference case; there is no reference equivalent. That is why
    `locus_copies == 1` can still produce a pair holding the reference — user visible as
    `ALLELE/Novel_gene_deletion`, whose phenotype then comes from the reference.
- **A missing token is two opposite things and only one of them is evidence.** A reference
  allele's defining variant absent from the pool means either that the alternate at that locus
  is homozygous, so there is no reference copy for it to sit on, or that the caller never
  looked. Both leave exactly the same gap. `locus_has_a_copy_number` (`phased.py:1422`) tells
  them apart by asking whether any token at that position carries a copy number, `NO_DATA`
  excepted — the line `variant_pool_numeric` already draws when it omits it.
  `no_defining_variant` (`:1471`) is gated on that and runs in both arms: nothing in it reads
  phase, and the `--phased` gate it used to carry was deciding something the flag is not about.
  - **The array is where reading absence as impossibility bites, and the default e2e set
    cannot see it.** Ungated, 547 array cells turned `X/X` into `Undetermined/Undetermined` —
    HPA3 413, GYPB 85, then HPA2, YT, HPA5, HPA1, LU, FY, A4GALT, KEL and JK — every one a
    locus the probe did not call. Against that, 281 pairs across the other inputs have a real
    call contradicting the reference. The two populations separate perfectly on the copy
    number question, and `1kg_microarray` is one of the two datasets `--full` exists for.
  - **`NO_COPIES` counts as evidence and `NO_DATA` does not.** Zero copies is a measurement of
    absence and a locus with no copies has no reference copy either; not measured is the
    absence of a measurement and says nothing. No input to hand has an instance of the
    `NO_COPIES` side, so a unit test is the only thing holding that choice in place.
  - **The locus test is a prefix match on `chrom:pos_`**, which is safe where the bare
    substring match elsewhere is not: `1:` cannot match `11:`, and `1:25408711_` cannot match
    `1:254087110_`. Do not simplify it to an `in`.
- **`FILTER` does not always mean call quality.** On some arrays PASS/FAIL is *probeset
  selection* — which of several probesets is the recommended one for a marker — so a FAIL row
  may be a perfectly good call. `only_keep_alleles_if_FILTER_PASS` will still drop it and revert
  to reference under a QC-sounding name. Since all QC is delegated upstream and `FILTER` is the
  only channel left, check what it means on a new input type before trusting it.
- **A `FILTER` value only costs anything if it lands on a defining variant of an allele
  that completes.** `only_keep_alleles_if_FILTER_PASS` walks the alleles that were *built*,
  and an allele is built only when every one of its defining variants is present — so a
  doubtful call that never completes an allele is never looked up. Measured over twelve
  sources and 9,281 samples: **91 values declared in headers, 17 present on a row at a
  database locus, 4 ever consulted, 1 (`LowQual`) ever excluding.** That funnel is why
  `filter_values.tsv` can be 26 rows against 91 declared values and still have no gap —
  and why counting unclassified values in a header says almost nothing about exposure.
  Two consequences worth knowing before touching the table:
  - **A value present but never consulted is a watch list entry, not a non-issue.** It is
    one sample away from mattering, and until it is classified it will exclude and revert
    to reference under a QC-sounding name.
  - **The table is keyed by value alone, and that is safe here but not by construction.**
    81 values are declared by more than one input form and none has a differing
    description. A value whose meaning differed between callers would need a second key
    and a loader that knows which caller wrote the file. Also: 9 of the 26 rows are not
    exercised by any input to hand, so the table's size is not evidence of coverage.
- **Bare `assert`s carry real validation** (`data_procesing.py:778, 1795`). They vanish under
  `python -O`. Prefer `BeyondLogicError` when replacing them.
- **An antithetical blood group is not always two subtypes, and the open question is narrower
  than it looks.** `antithetical_modifying_SNP_is_HOM` (`filters/geno.py`) once assumed one
  subtype per antigen — true for DO, FY, JK, KEL, KN and LU, **false for GYPB and A4GALT, which
  have three each**. GYPB carries `*03`, `*04` and `*06`.

  Three subtypes is not by itself a problem. The `len(putative_mod_SNPs) == 1` gate normally
  rules the case out, so the filter does nothing and the sample is fine. The guard therefore
  sits *after* the homozygosity test, where it only fires if a group with more than two
  subtypes also reaches the point of excluding pairs. **That combination has never been
  observed.** Five samples in a densely called short read cohort reach three subtypes and all
  pass, because the union was 7 or 9 rather than 1.

  So the live question is not "can a blood group have three subtypes" (yes) but "can three
  subtypes converge on a single modifying SNP that is homozygous" — which has **no real example
  to take to the lab scientists**, and is why the guard stays in place waiting for one rather
  than being widened or removed. If it ever fires, its `context` carries the sample, blood
  group, subtypes and SNP: that error is the whole evidence, so do not reduce it.

  Expect the population at risk to grow rather than shrink. `GYPB*06` needs five defining
  variants, so a probe based platform never builds it and the array and long read sets show
  nothing; ~2.4% of the short read cohort called it.
- **The database is hand curated, so whitespace turns up in it.** `prepare_db` calls
  `strip_stray_whitespace` on load for exactly this reason, and warns about what it trimmed so
  the fix lands in `db.tsv` rather than being papered over each run. It found 32 cells across
  `Coding`, `Protein`, `Sub_type` and `Note`. The one that mattered was `'GYPB*04 '` on seven
  rows, which forks one subtype into two everywhere `Sub_type` is used as a dictionary key.
  Also: since pandas 3 a text column reads back as `StringDtype`, **not** `object`, so
  `dtype == object` no longer finds string columns — that idiom silently does nothing now.
- **Test doubles drift silently.** `MockBG`, `MockBloodGroup` and `MockDb` in
  `tests/test_data_processing.py` are hand-rolled classes, not `spec=`'d mocks, so a new field
  on `BloodGroup` or `Db` breaks 28 tests at once with `AttributeError` rather than being
  inherited. Same class of problem as the local `Zygosity` stubs that shadow the real enum in
  three test files.
- **A multi-allelic row is one row but several pool entries, each with its own genotype.**
  `SmallVariantEncoder` emits one token per ALT, comma joined in ALT order, and
  `get_variants` fans them out and rewrites the genotype per alternate — `1/2` becomes `1/0`
  for the first token and `0/1` for the second, so `get_ref` and the phased filters only ever
  see `0` and `1`. Consequences: `vcf.variants` does **not** map one-to-one onto rows; a token
  recoding to hom ref is dropped, so an ALT the sample lacks is absent rather than present at
  zero; and the recode is applied *only* to rows with more than one ALT, so a single-ALT row
  with a GT like `2/2` still raises rather than being silently zeroed. Mixed SV/small
  multi-allelic rows (`ALT=<DEL>,G`) are **not** handled — one encoder claims the whole row.
- **A VCF can carry more than one row for one variant, and there are two unrelated reasons
  why.** Told apart by the FORMAT keys, which is the only reliable signal: two rows written
  by the same caller declare the same fields, two callers do not.
  - **Two callers, one variant.** A run emitting both a targeted caller's calls and the
    general caller's, with `FILTER` marking which row conflicts. Measured: 2,791 collisions
    in 591 of 967 samples, five RHCE tokens, one input form; the FORMAT keys differ in
    every one. `reconcile_duplicate_rows` (`vcf.py:618`) collapses these before anything
    reads them, keeping a phased genotype where either row has one and `;` joining the
    `FILTER` values.
  - **One caller, two events, one token.** `format_size` (`IO/encoders.py:128`) **floors**
    rather than rounds, so a token buckets every size within 1,000 bp — 479,720 collisions,
    up to **34 rows** on one token, 3,492 with different genotypes. **Deliberately not
    reconciled**: they are not two readings of one call, and `SvReader` takes its events off
    those rows. This is what the duplicate-genotype warning now reports.
    **Scoped and found harmless, which is worth knowing before anyone tries to fix it:** the
    token plays no part in matching — `select_best_per_vcf` groups by the *real*
    `(chrom, pos, end, svtype)` — and a matched allele's genotype comes from the matched
    event's own fields re-keyed under the database token (`main.py:483`), so nothing reads
    the collided entry. 3,188 lookups measured across every dataset with structural calls,
    none matching more than one row. Making the token precise is not a lever on anything.
  Reconciliation only touches single-base substitutions (`is_single_base_substitution`,
  `vcf.py:93`) — deliberately narrower than "not structural", because the structural
  reader's size threshold is `--min_size`, a runtime argument that layer never sees.
- **The `FILTER` lookup can match several rows, and it takes the first.**
  `filter_values_for` (`data_procesing.py:545`) is the one place it happens; it returns a
  list because the answer is genuinely plural, and the substring match is load-bearing for
  the multi-allelic fan-out, so exact equality would break it. Where the matched rows
  classify differently the verdict depends on file order, which
  `rows_disagree_about_exclusion` (`:592`) watches for. Zero on every input to hand — what
  holds it at zero is which column of `filter_values.tsv` a value sits in, not anything in
  the code.
- **A structural token's position is not a locus of the gene that carries the allele.** On a
  gene conversion allele the `DEL` is in this gene's coordinates and the paired `INS` is in
  the *donor's* — `RHD*01N.43` is `25272547_DEL_18244` in RHD plus `25402595_INS_18269` in
  RHCE, and `RHCE*03.02` is the exact mirror. So anything that groups positions by blood
  group from the raw tokens leaks across a paralogue pair. `get_loci_by_type` did, which put
  6 RHCE positions in RHD's set, 7 RHD positions in RHCE's, and made them share 8; since
  `locus_copies_for_bg` wants agreement across a gene, one borrowed position refused both
  genes on every sample. It now skips `_looks_like_sv_token`, which separates RHD
  (25272548-25328922) from RHCE (25362527-25420796) with nothing telling it where the genes
  are. Breakpoint coordinates are deliberately imprecise anyway — a GT landing on one is a
  coincidence, not evidence.
- **Nothing guarantees one VCF event per database SV definition.**
  `select_best_per_vcf` promises at most one match per *event*; the loop at
  `main.py:481-484` then writes `vcf.variants[f"{chrom}:{db.raw}"]` unconditionally, so two
  events best-matching one definition means **the last one silently wins the genotype**.
  This is in the path that decides the call, unlike the token collision above, and nothing
  reports it. Measured at 0 across 3,188 matched keys — but the shape that produces it is
  several near-identical events at one position, which is common on a merged or re-called
  file. Assume it becomes reachable if matching changes.
- **RH (RHD/RHCE) is long-read only** and explicitly unpolished. Short-read mismapping between
  the paralogs is considered intractable. Do not extend RH support to short read.
- Fuzzy SV matching was tuned on ~7 unique real SVs and is acknowledged as probably overfit.
  Treat its thresholds as provisional, not as validated constants.
- **In the alphanumeric phenotype file a `/` means two different things, and one of them
  is a delimiter.** `FUT1` and `FUT3` join a cell's alternatives with `/` *inside* a single
  entry (`choose_pheno.py:86`, `:149`), so `H+,Se+/H+,Se-` is two alternatives in one
  string where every other blood group would have written `H+,Se+ | H+,Se-`. Both touch
  `PhenoType.alphanumeric` only, so the numeric file never shows it. And the character
  cannot simply be split on, because `Ax/Aweak` and `U+alteredU/GPB` are single phenotype
  names that contain a literal `/`. So anything comparing phenotype cells as sets of
  ` | ` entries reads a FUT *narrowing* as a disagreement: on the long read phase pair
  that is 351 of the 381 the alphanumeric file appears to have, and 11 of the truth set's
  12, while the numeric file and the genotype file show neither. **Compare the genotype
  file.** Of the conflicts the two phenotype files hold once FUT is set aside, every one
  is already a conflict in the genotype file, so they add no signal of their own.
- **The three output files list a cell's alternatives in three different orders, and the
  lists are not even the same length.** The genotype cell is pair order,
  `",".join(bg.genotypes)` (`main.py:638`); both phenotype cells are
  `" | ".join(sorted(set(...)))` (`:680`, `:684`), which is alphabetical and
  **deduplicated**. So position *i* names a different pair in each file, and where two
  pairs share a phenotype string the phenotype list is *shorter* than the genotype list.
  **Per-cell lists are sets. Joining the three files by position is silently wrong**, and
  no promise of alignment can be added without dropping the deduplication, since four
  pairs really can all mean one phenotype. Measured over all nine datasets, of 19,913
  cells holding more than one genotype: **13,475** have a phenotype list of a different
  length, **9,193** have a genotype list in an order the sorted phenotype lists cannot
  match, **17,424** show at least one. `HG00125` JK is the clean case, six entries in
  every list, where position 0 gives one pair's genotype, that same pair's numeric
  phenotype and a *different* pair's alphanumeric phenotype. Two things that make it
  easy to miss: the `--debug` block prints all three lists one under the other
  (`record_data.py:213-215`), which is exactly the reading that does not hold; and the
  orders are **stable run to run**, so this is a systematic ordering difference rather
  than the `Pair.alleles` frozenset iteration hazard, which is separate and still latent.
  FUT1, FUT2 and FUT3 are excluded from those counts: the merge at `main.py:631-636`
  puts the other gene's pairs in the genotype cell, so the lists there are not describing
  the same thing to begin with.

## Design constraints inherited from history

The maintainer's note to developers: there was no upfront scope document; rules were teased
out of lab scientists iteratively, so the code and DB carry "evolutionary baggage". A
large-scale refactor was considered and rejected on ROI grounds — everything works and is
abstract enough to extend. The code is mix of functional and SOLID principles. The clean
code principles are also followed as much as possible without obsessing over it (ie best 
if a function doesn't have more than 3 parameters but there's a few with four). 

**Therefore: prefer additive changes over refactors.** Adding a parallel field or a new
resolution pass is welcome; restructuring core dataclasses or filters is too, but
only if strictly needed to avoid bad code/logic.

Also deliberately out of scope: QC metrics (use an upstream tool), novel variant reporting
(separate tool), hybrid alleles (need upstream assembly), Chido/Rodger (no in-house expertise).

## When changing filter logic

1. Reproduce with `--debug` and capture the log first.
2. Check the documentation PDF — many surprising behaviours are intended, and the docs
   describe them with worked examples. "Correct" is sometimes subjective here.
3. Add the case to the relevant `tests/test_*.py` before fixing.
4. If a filter's exclusion reason string changes, that string is user-visible and documented.
