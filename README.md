<table>
  <tr>
    <td>
      <h1>RBCeq2: blood group allele inference</h1>
    </td>
    <td align="right">
      <img src="images/Lifeblood-R_Primary_Keyline_RGB.jpg" alt="Lifeblood Logo" width="150">
    </td>
  </tr>
</table>

> [!WARNING]
> NOT FOR CLINICAL USE

## Version v2.4.4

RBCeq2 infers blood-group genotypes and phenotypes from supplied variant calls using
the curated ISBT allele definitions distributed with the package. Deterministic
rules select and exclude candidate alleles and pairs; `--debug` records the evidence
and named exclusions behind the result.

## Install and run

Python 3.12 or newer is required.

```bash
python -m pip install rbceq2==2.4.4
rbceq2 --version
rbceq2 --vcf sample.vcf.gz --out result --reference_genome GRCh38 --debug
```

Supply `GRCh37` or `GRCh38` to match the input coordinates. A single-sample VCF,
multi-sample VCF, or directory of separate `.vcf`/`.vcf.gz` files can be supplied.
Both ordinary gzip and uncompressed VCF input are read. For a directory:

```bash
rbceq2 --vcf samples/ --out cohort --reference_genome GRCh38 --processes 4 --debug
```

Use `rbceq2 --help` for all current options. Useful switches include:

| Option | Purpose |
|---|---|
| `--phased` | Use genotype phasing and phase-set information where available. |
| `--RH` | Include RHD/RHCE for the supported inputs described below. |
| `--HPAs` | Include human platelet antigen results. |
| `--no_filter` | Disable FILTER-based allele exclusion; inference checks still apply. |
| `--validate` | Run additional VCF format checks. Basic input checks always run. |
| `--PDFs` | Generate per-sample PDF reports. |
| `--min_size` | Minimum indel/SV size for fuzzy structural matching; default 10. |
| `--processes` | Worker processes; default 1. Large input loading also uses memory. |

## Input contract

Provide genotyped VCF calls from an upstream caller, with its calling and coverage
QC completed. **Native gVCFs are not supported in v2.4.4**: unspecified ALT
placeholders such as `<NON_REF>`/`<*>` and reference blocks require additional
handling. A gVCF can be read as text and still produce incomplete or misleading
inference; use the caller's final genotyped VCF.

GT must be the first FORMAT key when the row is used for inference. Called indices
must refer to REF (`0`) or an ALT declared by that row. Missing alleles use `.`.
Supported GT spellings use `/` or `|` between canonical integer indices; initial
phasing separators and leading-zero indices are explicitly unsupported.

Small-variant multiallelic records are processed into per-ALT calls internally.
Some reference-containing multiallelic configurations can remain partially
unresolved. Inspect `Undetermined` results and normalize such sites upstream where
appropriate. Mixed symbolic structural/small-variant records need caller-specific
assessment and are not a general supported input route.

Fully called genotypes naming more than two copies are interpreted only where the
per-ALT dosage is zero or all copies. Intermediate dosage receives the named
`get_ref/dosage_between_the_bounds` refusal. This is not general polyploid or pooled
sample inference.

### Chromosome copies and gene copies

These are different statements and lead to different output shapes:

- A valid called haploid GT outside PAR on X/Y supplies per-sample evidence for one
  chromosome in that region. No-calls supply no such evidence. PAR retains its
  two-chromosome interpretation.
- Gene copy number can be read from GT ploidy only when the upstream caller
  explicitly uses that convention. All called database small-variant positions
  reported for the gene must agree on haploidy; there may be only one reported
  position. Two chromosome slots remain, with a missing-gene-copy marker in one.
- GT shape alone cannot establish the caller's convention. Haploid assembly output,
  masked-PAR encodings, and other conventions must not be assumed equivalent to a
  gene-copy-number call. An uncalled GT is not evidence of zero gene copies.

Called structural deletions are evaluated against the curated database and can
constrain which chromosome carries an overlapping allele. The deletion's presence
and the gene-copy interpretation inferred from GT are kept distinct.

### RHD and RHCE

Use `--RH` for long-read VCFs with structural-variant information, or for callers
that explicitly encode RH gene copy number as GT ploidy. Ordinary short-read
variant calls without either source of information are unsupported for RH
inference because of the similarity between RHD and RHCE.

Compatibility depends on the caller's stated encoding convention. Synthetic
examples establish the software's behavior; they do not establish compatibility
with every caller or clinical validity.

### FILTER and missing data

By default, `PASS` and absent filtering (`.` or empty) do not exclude an allele.
Other FILTER values are classified: recognized values unrelated to call correctness
can be retained; call-correctness failures and unrecognized values can exclude.
Unrecognized values are named in warnings. `--no_filter` bypasses this FILTER-based
exclusion, while genotype validity and inference rules remain active.

A no-call such as `.`, `./.`, or a partial missing GT becomes `NO_DATA` and alleles
requiring that uncalled variant are excluded by name. All nonmissing indices must
still be valid for the row.

**A reference output does not establish sequencing coverage.** Missing records can
lead to reference defaults, and a no-call can still reach a reference fallback when
no candidate remains. Inspect the debug evidence and exclusions; a fallback is not
a measured homozygous-reference call.

## Outputs and interpretation

With `--out result`, the main outputs are:

- `result_geno.tsv`
- `result_pheno_numeric.tsv`
- `result_pheno_alphanumeric.tsv`
- `result_<run UUID>_log.txt`

The first TSV column identifies the sample; its heading contains the run UUID.
Cohort input uses the VCF sample names. Single-file and directory input uses the
filename with its final suffix removed, so `sample.vcf.gz` is named `sample.vcf`.

| Slot value | Meaning |
|---|---|
| A named allele | The selected database allele. |
| `-` | No second chromosome slot in the inferred single-copy region. |
| An absent-gene subtype, e.g. `RHAG*01N` | The database's designation for a missing gene copy without naming a specific breakpoint-defined allele. |
| `Novel_gene_deletion` | A missing gene copy for which no applicable database absence subtype is available. GT-based use depends on the caller convention above. |
| `Undetermined` | An allele slot that could not be named. Both slots can be undetermined when the blood group cannot be interpreted. |

`allele/Undetermined` retains one resolved slot. `Undetermined/Undetermined` can
represent a whole-group refusal; it does not establish the physical copy count.
Neither is a database allele.

Genotype alternatives are comma-separated. Phenotype alternatives are listed
independently and can be deduplicated, so entries must not be joined by position
between files. A slash can be part of a phenotype name.

Malformed retained GTs fail the affected sample and are reported with a named error;
the other samples in a cohort or directory can continue. Shared header/FORMAT
problems can prevent loading an input file. A nonzero process exit and the log must
be reviewed even when some result files were written.

## Database v2.5.1

The packaged [db.tsv](src/rbceq2/resources/db.tsv) is the source of allele definitions
used by inference. It is curated against ISBT definitions and includes explicit
modelling and nomenclature choices. Package and database versions are reported
separately by `--version`.

## Documentation and reporting issues

This README describes the v2.4.4 input and output contract. Longer PDF worked
examples attached to earlier releases are historical and may describe earlier
behavior; use the version-matched contract here when interpreting current output.

When reporting an issue, include the package/database versions, genome build,
command, and the complete affected sample/blood-group debug block where shareable.
State the caller and its genotype/copy-number convention. Expected answers can
require biological adjudication; a successful run or agreement with stored gold
does not establish clinical validity.
