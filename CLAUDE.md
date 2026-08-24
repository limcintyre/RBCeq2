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

`~/Dropbox/RBCeq2_related/scripts/e2e.py` runs the real CLI over nine datasets and diffs the
three output TSVs against gold standards in `~/Dropbox/rbceq2/e2e_gold/linux`:

| key | data | genome | flags |
|---|---|---|---|
| `1kg_microarray` | `~/Dropbox/vcfs/ALL_just_genes.vcf.gz` | GRCh37 | — |
| `ont_vienna_1kg_phased` | `~/Dropbox/vcfs/clair3_norm/` | GRCh38 | `--phased --RH` |
| `ont_vienna_1kg_unphased` | `~/Dropbox/vcfs/clair3_norm/` | GRCh38 | `--RH` |
| `public_truth_17_unphased` | `~/Dropbox/rbceq2/public_truth/combined_VCFs_uncompressed/` | GRCh38 | `--RH` |
| `public_truth_17_phased` | same | GRCh38 | `--phased --RH` |
| `dragen_per_sample` | `~/Dropbox/vcfs/DRAGEN/SNV_and_CNV_and_SV/per_sample/` | GRCh38 | `--RH` |
| `dragen_per_sample_phased` | same | GRCh38 | `--phased --RH` |
| `dragen_joint_967` | `.../cohort_967samples_SNV_SV_ROI50kb.vcf.gz` | GRCh38 | `--RH` |
| `dragen_joint_3209` | `~/Dropbox/vcfs/DRAGEN/SNV/cohort_3209samples_ROI50kb.vcf.gz` | GRCh38 | `--RH` |

```bash
python ~/Dropbox/RBCeq2_related/scripts/e2e.py --datasets public_truth_17_phased
```

**Seven of the nine run by default**, 4m 24s. `1kg_microarray` and `dragen_joint_3209`
are left out as the two slowest; `--full` runs all nine, 8m 40s, and either can still be
named directly. Every run adds `--HPAs --debug`.
Identical output to gold is a pass; otherwise `report_minimal_differences` prints
per-sample, per-column diffs (`e2e.py:407`), with comma-delimited fields compared as
unordered sets (`:310`).

`--filter-ab` runs every selected dataset a second time with filtering off and counts the
cells the tool names only in that arm (`e2e.py:618`, reported at `:956`). It doubles the run,
8m 17s against 4m 24s on the default seven, so it is off by default. See the bullet on checks
that need no gold below.

Datasets run two at a time with 8 workers each at `nice 10` (`--jobs`, `--processes`,
`--nice`). 16 workers on 16 cores is deliberate — a dataset reading a cohort VCF is
single threaded for a while and the other should use those cores. Use `--jobs 3` with
`--full`, where the third slot is worth 24 seconds; on the default seven it makes no
difference to wall time. `--jobs 1 --processes 12 --nice 0` is the old one-at-a-time
behaviour. Concurrency does not change output: every arm of the sweep produced all 27
TSVs byte identical to the serial run.

**A run ends with a summary of what differs from gold**, so answering "did anything
change" does not mean scrolling back through tens of thousands of lines. It is a tally
of differences with a pointer to them, never a verdict — see the next bullet.

How it is used, so agents read the results correctly:

- **The maintainer reviews every gold-vs-new discrepancy by hand.** The script is a difference
  *reporter*, not a pass/fail gate — the exit code carries no signal, and the end-of-run
  summary counts differences rather than declaring an outcome. Never write "e2e passed"; e2e
  produces a diff, and a human adjudicates it. A difference may well be an improvement, and an
  agreement may be two copies of the same defect — see "e2e cannot see a defect that was
  present when gold was made" in `TODO.md`.
- **It runs the installed `rbceq2` console script** (`e2e.py:244-245`), not the working tree, so
  the active env's install has to be current for a change to show up in the output. Check this
  before reading any e2e result, and before making gold from one — gold built from a stale
  install bakes in whatever that install was. This is not hypothetical: on 2026-08-24 the
  install was found two commits behind, and the day's gold had been built from it. One command
  settles it, and it is worth running rather than assuming:
  `diff -rq src/rbceq2 "$(python -c 'import rbceq2,os;print(os.path.dirname(rbceq2.__file__))')"`
- **Gold standards are platform- and version-specific** (`.../e2e_gold/linux`, a given DB
  version). When a change is *supposed* to alter output, the gold files need regenerating —
  that is the maintainer's call, never an agent's.
- **A key with no gold reports that and moves on** (`e2e.py:551`). It is a normal state, not a
  failure: the run still produces the output a gold would be made from. All nine have gold as
  of 2026-08-20; this exists so that adding a tenth does not end the run at the new key.
- **The default set has no GRCh37 and no array.** `1kg_microarray` is the only dataset that is
  either, and it is one of the two the default leaves out. The `FILTER` landmine below —
  PASS/FAIL meaning probeset selection rather than call quality — is about arrays specifically,
  and nothing in the default set can see it. Use `--full` before trusting a change that touches
  build handling, `FILTER`, or lane variants.
- **The last four keys carry a check that needs no gold.** `dragen_per_sample` and
  `dragen_joint_967` are the same 967 people, and `dragen_per_sample_phased` is the first of
  those read a second way. Joint calling changes the encoding, not the biology, and phase may
  narrow a call but must never change one — so a cell where neither arm's set contains the
  other is the tool contradicting itself, which is a stronger signal than a gold diff because
  it needs no gold to interpret. Neither comparison is automated; both are read off the three
  TSVs by hand, comma fields as unordered sets. Last measured over 85,096 cells: the phase
  comparison is 84,196 agree / 900 narrowed / **0 conflict**, and a conflict there is a defect.
  The encoding comparison is 85,002 / 60 / **34**, and that 34 is not the tool's fault — the
  two files disagree about the genotype at some sites, all 34 are RHCE, so read a *change* in
  the number rather than the number itself.
- **`--filter-ab` is a third check that needs no gold, and the only automated one.** It runs
  each dataset again with filtering off and counts cells the tool declines to name in the
  first arm and names in the second. The flag has one behavioural site in the tool
  (`main.py:493`), so a difference is attributable to that filter and nothing else, and both
  arms are the same input read by the same code — which is why it can see a defect that was
  already in the output when gold was made from it.

  **A non-zero count is not a defect.** Filtering removing a call is what filtering is for.
  Measured over all nine datasets on 2026-08-24, after `f33f1c3`: **no whole cells and nine
  single slots**, in four samples and two blood groups, all long read —
  `HG01872`/`NA18544` GYPA and `HG03730`/`HG03886` RHCE, each in both long-read arms, plus
  `NA18571` RHCE in the phased arm only. Every one is the tool behaving as designed; they are
  the two `cant_name_second_slot_cuz_*` filters doing their job, and each filter's docstring
  example is one of them. So read a *change* in the set, not the count: a cell arriving or
  leaving is the signal. Traces for all of them are in the working directory.

  The whole-cell column was 1 before `f33f1c3` and is the check's one scalp so far: it found
  `NA18571`, which is now `RHCE*01/Undetermined` rather than `Undetermined/Undetermined`.

  In every one of the five, the excluding value is `LowQual` and nothing else, on a call the
  caller gave GQ 0 or 1 at 3-11x depth. The `FILTER` decisions themselves look right; what is
  worth knowing is that the PASS calls beside them are GQ 2-7, so on this input the two sides
  of the line are closer together than the names suggest.

  **The nine are recorded in `FILTER_AB_EXPECTED` and the run says what moved**, so the count
  is not a number anyone has to remember. It records membership — dataset, sample, blood group,
  and which of the three ways the arms disagreed — and deliberately not the genotypes: a cell
  whose *value* changes is a gold difference and the gold diff reports it. The two checks
  divide the work, gold owning the values and this owning the membership, which is the half
  gold is blind to. Regenerate it the way gold is regenerated: when a change is supposed to
  move it, and by the maintainer.

  Two things it deliberately does not do. It reads only the genotype TSV, because a refusal
  propagates into both phenotype files and counting those would count one cell three times.
  And it only reports cells the unfiltered arm names *in full*, so where both arms are partial
  there is no signal — a floor, not a zero, and a real one: of the 23 cells `f33f1c3` fixed,
  22 were invisible here because both arms refused them, and only `NA18571` ever showed up. It
  also reports the reverse direction separately, which is not a defect: the unfiltered pool
  holds variants the caller doubted, and one of those can contradict a reference the filtered
  pool left standing.

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

- **A haploid GT has three readings, and two of them are live.** The same `1` means different
  things in different places and the results have **different shapes**, so `get_ref`
  (`data_procesing.py:1363`) decides by which of two counts is 1 rather than by the GT:

  | which count is 1 | what it means | allele slots | second slot |
  |---|---|---|---|
  | `chrom_copies` | outside PAR on X/Y; the sample has one chromosome | **1** | `-` (`HAPLOID_SECOND_SLOT`) |
  | `locus_copies` | the caller encoded gene copy number as GT ploidy | **2** | the DB's subtype for a missing copy, else `?` (`UNNAMED_SECOND_SLOT`) |
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
  chromosome, `?` says there is one and it carries no copy of the gene, `Undetermined` says it
  carries a copy whose allele the database cannot name. Collapsing any pair of them undoes the
  distinction `chrom_copies`/`locus_copies` exists to draw, and pairing with the reference
  instead — the option that looks tidiest — asserts wildtype on a chromosome there is positive
  evidence against. All three are user visible and documented; see `constants.py:584-615`.
- **`./.` is overloaded two ways**: a genuine no-call and, from copy-number-aware callers, zero
  copies. These mean opposite things. The synthesised lane row is *not* one of them any more —
  it carries `SYNTHESISED_HOM_REF_GT` (`constants.py:23`), which is the only legitimate
  assertion of wildtype in the codebase.
- **`Zygosity` is a projection of two numbers, not one label.** `HEM` is one copy of the locus
  carrying one copy of the token; `NO_COPIES` is zero of both. Neither says how many allele
  slots the *result* has — that is `BloodGroup.chrom_copies`, a third number, and it is never
  changed by a deletion. Reading pair shape off a per-token zygosity is the bug issue #40 was.
- **`FILTER` does not always mean call quality.** On some arrays PASS/FAIL is *probeset
  selection* — which of several probesets is the recommended one for a marker — so a FAIL row
  may be a perfectly good call. `only_keep_alleles_if_FILTER_PASS` will still drop it and revert
  to reference under a QC-sounding name. Since all QC is delegated upstream and `FILTER` is the
  only channel left, check what it means on a new input type before trusting it.
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
