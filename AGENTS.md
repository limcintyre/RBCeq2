# AGENTS.md

Standing instructions for RBCeq2. Read once before work. `CLAUDE.md` is a pointer
here, not a second policy. When resuming, read `NEXT_SESSION.md` in the session
working directory, then only the task's named evidence and relevant reference
sections below. Do not load every linked document or reread unchanged instructions.

## Project and hard rules

RBCeq2 (`rbceq2`, Python >=3.12) reads VCFs and infers blood group genotypes and
phenotypes against the ISBT allele database curated by Australian Red Cross
Lifeblood. Its core promise is deterministic inference with an auditable exclusion
trail.

- **NOT FOR CLINICAL USE.** Preserve the banner in the README, documentation, and
  runtime output. Do not soften it or imply clinical validity.
- Same VCF and database must give the same result. No probabilistic filtering,
  Bayesian inference, likelihood scoring, ML, or randomness.
- Every removed allele or pair must be recorded with its named exclusion reason.
  Silent drops are bugs; filter names and reason strings are user-visible behavior.
- The curated `src/rbceq2/resources/db.tsv` is the source of truth for biology.
  Investigate code/database disagreements there; do not invent allele definitions
  in code. Correct erroneous database definitions in the database.
- Prefer explicit, unambiguous allele definitions, phenotype strings, and logs.
  Use Google-style docstrings for new prose; preserve existing real log examples
  even when their formatting is imperfect. Keep notes and comments technology-agnostic.
- **Version stays 2.4.4 until the maintainer explicitly agrees otherwise.** Do not
  bump the package or database version, including implied future releases in
  docstrings. Date a behavior change by what changed. If a package version change
  is agreed, check all four version locations; database versioning is separate.

## Working agreement: propose, never apply

- Confirm repository path, session working directory, and HEAD before work. The
  repository is read-only for agents; the maintainer may edit it during a run.
- Give the whole plan as dot points first, including decisions the maintainer must
  make. Propose one small step at a time; the maintainer applies it before the next.
- Write complete proposed files in the session working directory with the same
  basename. Hand over `code --diff <repo path> <working dir path>`, repository first.
  For new files, identify the destination and provide a command for the maintainer
  to copy them. Anchor findings to `file:line` and keep changes self-contained.
- All findings, scripts, logs, snapshots, and proposals go in the session working
  directory or `/tmp`. **Never delete anything and never write elsewhere.** The
  repository, golds under `~/Dropbox`, and `e2e.py` changes come to the maintainer as
  commands. Prefer a `for` loop for repeated operations. If asked to cross this
  rule, say so first and check the command on a copy or with `echo` before handing
  it over. Account for test cleanup and install side effects when choosing commands.
- Do not overwrite the maintainer's edits or clean their worktree. Preserve the
  previous working proposal before replacing it. Do not regenerate golds; the
  maintainer applies diffs, reviews output changes, and regenerates golds.
- End each step with a concise git commit message: subject plus **one paragraph**.
  Name only files the commit contains; session notes such as `TODO.md` and
  `NEXT_SESSION.md` are not repository files.

## Delegate bounded work

Use subagents for concrete independent source exploration, implementation in
proposal copies, test execution, and log analysis. Keep requirements, decisions,
integration, and final review in the main conversation.

- Give narrow task briefs, relevant paths, and all applicable constraints. Avoid
  copying the whole conversation or asking agents to reread duplicate policies.
- Assign one owner per proposed file. Request concise findings with `file:line`
  and evidence paths; save full logs in the working directory rather than returning
  them to the main conversation.
- Check CPU load before launching jobs and coordinate every agent's jobs as one
  workload. This is the maintainer's workstation; they also run e2e here.
- All subagents follow the same proposal, location, deletion, and version rules.
  Delegation grants no permission to edit the repository or other restricted paths.

## Validation and evidence

- Reproduce a bug with `--debug` before proposing a source correction. For filter
  changes, read the relevant documentation PDF examples and check/add the relevant
  unittest case before fixing. Reuse already-applied regression cases rather than
  proposing them again. Check related pitfalls in `RBC_EQ_REFERENCE.md` before
  touching that code.
- Supply the **full raw debug block** from the sample/BG line through the `______`
  terminator. No excerpts or reformatting; timestamps stay intact except that the
  leading timestamp may be removed. Save it in the working directory and tell the
  maintainer the path. Tool output alone does not reach them.
- For proposal before/after comparisons, build baselines from a fixed commit or
  preserved snapshot, never the live working tree. Record source/input provenance
  and compare against that current-code baseline, not gold. Compare comma-delimited
  TSV fields as unordered sets; phenotype `/` is not a general delimiter.
- Use `unittest`, not pytest. The suite is a smoke check, not acceptance. State
  which tests ran and which end-to-end measurements or maintainer reviews remain.
- The maintainer runs e2e. Its gold differences need human adjudication; exit code
  and summary counts are not verdicts. Never write "e2e passed". Before interpreting
  e2e results, verify which installed `rbceq2` code ran. Changes to build handling,
  `FILTER`, or lane variants require coverage of the GRCh37 array omitted by the
  default seven datasets; see `RBC_EQ_TESTING.md` for `--full` and comparison rules.
- The cohort whose caller encodes gene copy number as GT ploidy is private and not
  shareable. Local VCFs carry haploid GTs only on non-PAR X. Where that unavailable
  input is needed, reason from the database and unit/synthetic tests and say so.

## Code and scope

Preserve intentional filter names and abbreviations (`cuz`, `SNP`) and the typo in
`core_logic/data_procesing.py`. Pipeline order is load-bearing; use
`@apply_to_dict_values` for new stages. Log with `loguru`; do not downgrade warnings
to debug to quiet them. Prefer additive changes; restructure core dataclasses or
filters only when strictly needed for sound code or logic.

Keep chromosome copies, locus copies, and token copies distinct. Preserve the
meaning of `-`, `Novel_gene_deletion`, and `Undetermined`. Do not change the reference
exemption in `pair_can_exist`; read its detailed rationale before related work.
RH handling remains long-read only. QC metrics, novel variant reporting, hybrid
alleles, and Chido/Rodger are deliberately outside scope.

## Read references only when relevant

Search the topic index, then read the matching sections. Historical measurements,
line numbers, and implementation status must be checked against current code.

| Work | Reference |
|---|---|
| Commands, e2e datasets, installed-code checks, gold interpretation, phase/encoding comparisons, filter A/B | [RBC_EQ_TESTING.md](RBC_EQ_TESTING.md) |
| Repository map, pipeline, glossary, naming and logging conventions | [RBC_EQ_REFERENCE.md](RBC_EQ_REFERENCE.md) |
| Ploidy, copy counts, reference pairing, absent tokens, FILTER, database parsing, VCF rows, structural variants, phenotype output | Relevant topics under Known landmines in [RBC_EQ_REFERENCE.md](RBC_EQ_REFERENCE.md) |
| Current step, accepted results, pending decisions and evidence | `NEXT_SESSION.md` in the session working directory |

Common commands, from an appropriate snapshot or by the maintainer as required by
write restrictions:

```bash
python -m unittest discover tests
coverage run -m unittest discover tests && coverage report
flake8 src/
rbceq2 --vcf x.vcf.gz --out y --reference_genome GRCh38 --debug
```
