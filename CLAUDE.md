# CLAUDE.md

Context for AI agents working in this repo. Read before making changes.

## What this is

RBCeq2 reads genomic variants (VCF) and infers blood group (BG) genotypes and phenotypes
against ISBT-defined alleles. It builds every possible allele from the observed variants, then
filters candidates with explicit logic checks until only viable pairs remain. The value
proposition is the **auditable trail** — every exclusion is recorded with a named reason.

Package `rbceq2` v2.4.3 · database v2.5.0 · Python ≥3.12 · maintained by Australian Red Cross
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
6. **Be explicit over terse.** Use google style doc strings, however, many of the examples 
   in the existing functions have copy paste from logs and this is not always formatted perfectly,
   Do no change these - they are real examples of the function (normally a filter).

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

1. **Read** — VCF parsed; SVs extracted (`SvReader`) and fuzzy-matched to DB tokens (`SvMatcher`)
2. **Raw alleles** — every DB allele whose defining variants are all present
3. **Variant pool** — `make_variant_pool` assigns zygosity per variant token
4. **Large-indel adjustment** — `modify_variant_pool_if_large_indel` rewrites HOM→HEM inside deletions
5. **Phasing** — phase pool built from raw GT strings; unphased alleles removed if `--phased`
6. **Pairs** — all viable pairs generated from the pool
7. **Filters** — ~23 named filters, plus Knops co-existing logic
8. **Phenotype** — antigen merge, modifier weighting, antithetical reconciliation
9. **Output** — 3 TSVs (genotype, phenotype numeric, phenotype alphanumeric), log, optional PDFs

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

- **`get_ref` (`data_procesing.py:856`) assumes a 3-character GT** and never splits on `/` or
  `|`. Haploid GTs crash it. It also does `.replace(".", "0")`, which converts absence of data
  into a positive assertion of wildtype. Open issue #40.
- **`./.` is overloaded three ways**: genuine no-call, the synthesised `HOM_REF_DUMMY_QUAL`
  lane row, and (from copy-number-aware callers) zero copies. These mean opposite things.
- **`remove_home_ref` (`vcf.py:120`) matches `"0/0"` only** — misses `0|0` and haploid `0`.
- **`filter_vcf_metrics`** evaluates `float(variant_metrics[variant][metric_name])`
  unconditionally *before* the `if microarray:` branch, so the microarray path still raises on
  missing DP/GQ. `--microarray` is currently a one-line depth stub carrying a `TODO !!`.
- **`Zygosity.REF`** is defined and given a `len_dict` value but is never produced or compared
  anywhere. Probably a vestigial "zero copies" cell.
- **Bare `assert`s carry real validation** (`data_procesing.py:535, 549, 878`). They vanish
  under `python -O`. Prefer `BeyondLogicError` when replacing them.
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