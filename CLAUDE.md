# CLAUDE.md

Context for AI agents working in this repo. Read before making changes.

## What this is

RBCeq2 reads genomic variants (VCF) and infers blood group (BG) genotypes and phenotypes
against ISBT-defined alleles. It builds every possible allele from the observed variants, then
filters candidates with explicit logic checks until only viable pairs remain. The value
proposition is the **auditable trail** — every exclusion is recorded with a named reason.

Package `rbceq2` v2.4.4 · database v2.5.0 · Python ≥3.12 · maintained by Australian Red Cross
Lifeblood.

## Hard rules

1. **NOT FOR CLINICAL USE.** This banner appears in the README, the docs, and at runtime. Do
   not remove it, soften it, or add language implying clinical validity.
2. **Deterministic, not probabilistic.** Same VCF + same DB ⇒ same result, always. No Bayesian
   filtering, no likelihood scoring, no ML, no randomness. This is a deliberate design choice
   ("not an AI black box"); a PR that adds probabilistic filtering breaks the core promise.
3. **Every exclusion must be recorded.** When a filter removes an allele or a pair, it goes
   into the exclusion log with the filter's name. Silent drops are bugs.
4. **The DB is the source of truth for biology.** `src/rbceq2/resources/db.tsv` is manually
   curated and is the origin of most errors. Do not invent allele definitions in code. If code
   and DB disagree, fix the DB.
5. **Be explicit over terse.** Verbose unambiguous output is preferred to concise ambiguous
   output. This applies to allele definitions, phenotype strings, and log messages.
6. **Real examples in doc strs.** Use google style doc strings, however, many of the examples 
   in the existing functions have copy paste from logs and this is not always formatted perfectly,
   Do no change these - they are real examples of the function (normally a filter).

## Working agreement

Agents **propose; the maintainer applies.** This repo is the maintainer's primary work.

1. **Do not edit files in this repo.** Output suggested changes as code blocks in chat, or as
   files written to the session working directory, for manual review and paste. This is
   deliberate: it keeps a human on every line that lands.
2. **Write all output to the working directory.** Findings docs, notes, analyses, throwaway
   scripts, and proposed versions of repo files all go there — never into the repo tree.
   `git status` should stay clean unless the maintainer dirtied it.
3. **Anchor every claim to `file:line`** so a suggestion can be found by hand. Prefer small
   self-contained blocks over large rewrites — they are being pasted, not applied.

## Commands

```bash
pip install -e .                      # dev install
python -m unittest discover tests     # tests (unittest, NOT pytest)
coverage run -m unittest discover tests && coverage report
flake8 src/
rbceq2 --vcf x.vcf.gz --out y --reference_genome GRCh38 --debug
```

No CI workflows exist. There is no Makefile. `--debug` produces the full per-sample filter
trace and is the primary debugging tool — reproduce a bug with `--debug` before changing code.

**The `unittest` suite is a smoke check, not the acceptance gate.** The end-to-end tests are
the real gate and the maintainer runs them manually after any change. Never report work as
verified on the strength of `python -m unittest discover tests` alone — state which unit tests
ran and that e2e is still outstanding.

### End-to-end tests

`~/Dropbox/RBCeq2_related/scripts/e2e.py` runs the real CLI over five datasets and diffs the
three output TSVs against gold standards in `~/Dropbox/rbceq2/e2e_gold/linux`:

| key | data | genome | flags |
|---|---|---|---|
| `1kg_microarray` | `~/Dropbox/vcfs/ALL_just_genes.vcf.gz` | GRCh37 | — |
| `ont_vienna_1kg_phased` | `~/Dropbox/vcfs/clair3_norm/` | GRCh38 | `--phased --RH` |
| `ont_vienna_1kg_unphased` | `~/Dropbox/vcfs/clair3_norm/` | GRCh38 | `--RH` |
| `public_truth_17_unphased` | `~/Dropbox/rbceq2/public_truth/combined_VCFs_uncompressed/` | GRCh38 | `--RH` |
| `public_truth_17_phased` | same | GRCh38 | `--phased --RH` |

```bash
python ~/Dropbox/RBCeq2_related/scripts/e2e.py --datasets public_truth_17_phased
```

All five run by default. Every run adds `--HPAs --debug --processes 12`. Identical output to
gold is a pass; otherwise `report_minimal_differences` prints per-sample, per-column diffs
(`e2e.py:163`), with comma-delimited fields compared as unordered sets (`:115`).

How it is used, so agents read the results correctly:

- **The maintainer reviews every gold-vs-new discrepancy by hand.** The script is a difference
  *reporter*, not a pass/fail gate — `validate_outputs` always returns `True` (`e2e.py:264`),
  so the exit code carries no signal. Never write "e2e passed"; e2e produces a diff, and a
  human adjudicates it.
- **It runs the installed `rbceq2` console script** (`e2e.py:60-61`), not the working tree, so
  the active env's install has to be current for a change to show up in the output.
- **Gold standards are platform- and version-specific** (`.../e2e_gold/linux`, a given DB
  version). When a change is *supposed* to alter output, the gold files need regenerating —
  that is the maintainer's call, never an agent's.

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
    vcf.py                   VCF/gVCF reading, phase sets, lane variants
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

- **A haploid GT is only a ploidy statement outside PAR on X/Y.** Anywhere else — inside PAR, or
  on an autosome — it is a statement about *locus* copy number and has a different correct
  answer, so `get_ref` (`data_procesing.py:1085`) rejects it by name rather than guessing. Some
  platforms encode gene copy number as GT ploidy, so a haploid GT across a gene with a
  whole-gene deletion allele means one copy of *that gene*, which is still two allele slots —
  one of them the deletion allele — **not** one slot. Mapping that needs
  `ploidy_state_table.md` decisions 2 and 4, which are unmade. Do not route it through the
  haploid path.
- **`./.` is overloaded two ways**: a genuine no-call and, from copy-number-aware callers, zero
  copies. These mean opposite things. The synthesised lane row is *not* one of them any more —
  it carries `SYNTHESISED_HOM_REF_GT` (`constants.py:23`), which is the only legitimate
  assertion of wildtype in the codebase.
- **`Zygosity` is a projection of two numbers, not one label.** `HEM` is one copy of the locus
  carrying one copy of the token; `NO_COPIES` is zero of both. Neither says how many allele
  slots the *result* has — that is `BloodGroup.chrom_copies`, a third number, and it is never
  changed by a deletion. Reading pair shape off a per-token zygosity is the bug issue #40 was.
- **`filter_vcf_metrics`, `remove_alleles_with_low_read_depth` and
  `remove_alleles_with_low_base_quality` are dead code** — defined, documented, tested, and
  never called from anywhere in `src/`. There is no `--microarray`, `--min_read_depth` or
  `--min_base_quality` flag; all QC is deliberately delegated upstream and enforced only through
  `FILTER == PASS`. The tests are the sole reason they still pass. Do not wire them back in
  without asking.
- **`FILTER` does not always mean call quality.** On some arrays PASS/FAIL is *probeset
  selection* — which of several probesets is the recommended one for a marker — so a FAIL row
  may be a perfectly good call. `only_keep_alleles_if_FILTER_PASS` will still drop it and revert
  to reference under a QC-sounding name. Since all QC is delegated upstream and `FILTER` is the
  only channel left, check what it means on a new input type before trusting it.
- **Bare `assert`s carry real validation** (`data_procesing.py:715, 1792`). They vanish under
  `python -O`. Prefer `BeyondLogicError` when replacing them.
- **Test doubles drift silently.** `MockBG`, `MockBloodGroup` and `MockDb` in
  `tests/test_data_processing.py` are hand-rolled classes, not `spec=`'d mocks, so a new field
  on `BloodGroup` or `Db` breaks 28 tests at once with `AttributeError` rather than being
  inherited. Same class of problem as the local `Zygosity` stubs that shadow the real enum in
  three test files.
- **RH (RHD/RHCE) is long-read only** and explicitly unpolished. Short-read mismapping between
  the paralogs is considered intractable. Do not extend RH support to short read.
- Fuzzy SV matching was tuned on ~7 unique real SVs and is acknowledged as probably overfit.
  Treat its thresholds as provisional, not as validated constants.

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
