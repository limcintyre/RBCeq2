# RBCeq2 testing reference

The body below was moved verbatim from the AGENTS.md snapshot saved on 2026-09-07.
Read the relevant sections on demand; current operating rules are authoritative in AGENTS.md.
Measurements, timings, commit/version mentions, and source line numbers record their original
observations. Check current code and data before treating those statements as current.
For e2e reporting, follow AGENTS.md: report differences for human review; never write "e2e passed".

## Topic index

| Read for | Line in this file |
|---|---|
| Commands; debug reproduction; unit-test limits | 22 |
| E2e datasets, command, flags, and workload | 40 |
| Gold differences, installed-code check, and gold ownership | 87 |
| GRCh37, arrays, and when to use `--full` | 106 |
| Phase/encoding comparisons, covering semantics, and historical cases | 111 |
| `--filter-ab`, expected membership, and limits | 183 |

<!-- Original AGENTS.md section content begins below; retained byte for byte. -->

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
- **Seven of the nine keys carry a check that needs no gold, and it is scripted.** Every
  phased key has an unphased twin, and `dragen_per_sample` and `dragen_joint_967` are the
  same 967 people read two ways. Joint calling changes the encoding, not the biology, and
  phase may narrow a call but must never change one — so a cell where neither arm's answer
  covers the other is worth looking at, and it needs no gold to interpret, which makes it a
  stronger signal than a gold diff. `survey_phase_comparison.py` runs all four comparisons
  off the genotype TSVs, comma fields as unordered sets; it used to be done by hand and it
  reproduced every number that had been derived that way.

  **Genotype file only.** A refusal propagates into both phenotype files, so counting those
  would count one cell three times, and the alphanumeric file's `/` is not a delimiter.

  Read under the **covering** reading: `Undetermined` is not an allele name, it is the tool
  declining to name that slot, so `X/Undetermined` covers `X/anything` and the phased arm has
  refused rather than disagreed. `-` and `Novel_gene_deletion` deliberately do **not** subsume — one
  says there
  is no second chromosome and the other that there is one carrying no copy of the gene, and a
  named allele contradicts those rather than refining them, so letting them cover would
  launder real conflicts.

  | comparison | cells | agree | narrowed | refused | **conflict** |
  |---|---|---|---|---|---|
  | short read phase | 85,096 | 84,951 | 145 | 0 | **0** |
  | long read phase | 57,904 | 55,343 | 2,530 | 2 | **29** |
  | truth set phase | 1,496 | 1,412 | 83 | 0 | **1** |
  | encoding | 85,096 | 85,002 | 20 | 1 | **73** |

  **A conflict is a defect on the phase pairs**, and the short read pair's 0 is the one that
  has always been held there. It was 84,196 / 900 / 0 before the three filters that read copy
  number where the phased arm reads phase — the two `..._shared_variant_has_too_few_copies`
  and `no_defining_variant` ungated — which between them close 755 of those narrowings by
  reaching the same verdict without phase.

  **The long read 29 and the truth set 1 are adjudicated and parked**, so read a *change*
  rather than the number. They are all one thing: 27 RHCE plus the truth set's 1, and 2 RAPH,
  where the unphased arm reaches an answer phase rules out — see below. Not a defect in either
  arm.

  **A co-existing slot is compared as a set of alleles, not as a string**, and that is a fix
  to the check rather than a reading of it. Knops writes more than one allele on one
  chromosome joined with `+`, and phase can only ever *remove* one from that group, so a
  phased slot holding fewer alleles is a narrowing. Compared as strings, `KN*01.07` against
  `KN*01.07+KN*01.12` looked like a different allele: that was 13 of the long read conflicts,
  now 12 narrowings and one refusal, and `NA18499` is the worked case. The direction is
  asymmetric on purpose — a phased slot carrying an allele the unphased slot never had is
  still a conflict, because phase *adding* an allele to a chromosome is the thing worth
  catching. Off Knops every slot holds one allele and this is equality, as it always was.

  **The refused column is not a quiet bucket.** Ten cells sat in it that were a real defect:
  `ref_not_phased` was discarding a chromosome phase had settled, fixed by
  `cant_name_second_slot_cuz_ref_not_phased`. A cell arriving there is still a prompt to look.

  **The invariant assumes the unphased arm's candidates contain the truth, and on RHCE they
  do not.** Two of the four filters that write a single named slot run unphased, so both arms
  can say `X/Undetermined` — but all four fire only where *nothing* paired. On these cells the
  unphased arm reaches a pairing that balances the pool, so it never gets there, while phase
  rules that pairing out and leaves the tool naming one chromosome and refusing the other.
  `HG00365` is the shape: unphased `RHCE*01.01/RHCE*01.36`, phased `RHCE*01/Undetermined`,
  and phase puts one token of each unphased allele on each chromosome. So the unphased answer
  is not less specific, it is wrong, and no reading of the comparison can make the two agree.
  That is what the 27 RHCE are, and it is why the fix above *raised* the long read count
  while improving the output. The 2 RAPH are the same shape on a compound heterozygote:
  `RAPH*01.-01.01` needs two variants and phase puts them in trans, so neither chromosome is
  a named variant allele and `RAPH*01/RAPH*01`, MER2+, is right — the unphased arm has to
  assume cis and gets it wrong.

  The encoding comparison's 73 is not the tool's fault either — the two files disagree about
  the genotype at some sites, all of them RHCE. Those same two filters took it from 60 / 34,
  and the 41 it gained are one shape: the per-sample file calls `1:25408711 G>A` where the
  joint file calls the same sample `0/0`, so a pair impossible in one arm is genuinely possible
  in the other. Shrinking the per-sample answer broke a containment that had been accidental —
  the disagreement was already there, and this is the case the sentence above is warning about.
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

