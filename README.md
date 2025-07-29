<table>
  <tr>
    <td>
      <h1>RBCeq2 documentation</h1>
    </td>
    <td align="right">
      <img src="images/Lifeblood-R_Primary_Keyline_RGB.jpg" alt="Lifeblood Logo" width="150">
    </td>
  </tr>
</table>


<div style="color: red;">

> **Warning**
> NOT FOR CLINICAL USE

</div>

RBCeq2 reads in genomic variant data in the form of variant call files (VCF) and outputs blood group (BG) genotypes and phenotypes. The user docs (Word) explains how RBCeq2 constructs possible allele combinations and then filters them until only possible genotype/phenotype combinations remain. Internal to the tool, all variants are based on human genome release GRCh37 or GRCh38, and in the context of a variant ‘_ref’ means no change to the reference nucleotide for the associated position.

## Overview

At the highest level RBCeq2 finds all possible alleles, then filters out those that fail certain logic checks. This allows for an auditable trail of why it has reached a certain result. Every effort has been made to be explicit both in encoding alleles in our database and while writing code. This results in verbose but unambiguous results. Some tools in this space employ Baysian or other likihood methods to filter results – we chose not to do that initially but might add it as an optional feature in the future. Last, some liberties have been taken to standardise syntax and nomenclature across blood groups.

This initial release of RBCeq2 is focused on perfecting the calling of International Society for Blood Transfusion (ISBT) defined BG alleles from simple variants (i.e single nucleotide variant -SNVs/small INDEL) that can be found in standard short read or microarray derived VCFs. Further, it supports the use of long read VCFs (i.e large indels and phased data). However, these features are not as polished.

## Bugs

This software is extensively tested and accurately reports genotypes/phenotypes based on our inhouse definitions of the ‘correct’ answer, however, there are some examples where the ‘correct’ answer is subjective. These docs are detailed – if you find what you think is a bug in the results from RBCeq2 please take the time to understand if it inline with what we intended or not. We will endeavor to fix any black and white bugs in less than one week. Most of these will be rare variants that are encoded wrong in our database. Further, we value any and all feedback and feature requests.

## How To

Install via pip (python3.12+) or clone the git repository:

```bash
pip install RBCeq2

rbceq2 -h

usage: py main.py --vcf example_multi_sample.vcf.gz --out example --reference_genome GRCh37

options:
  -h, --help            show this help message and exit
  -v, --version         Show program's version number and exit.
  --vcf VCF             Path to VCF file/s. Give a folder if you want to pass multiple separate files (file names must end in .vcf or .vcf.gz), or alternatively give a file if using a multi-sample VCF. (default: None)
  --out OUT             Prefix for output files (default: None)
  --depth DEPTH         Minimum number of reads for a variant (default: 10)
  --quality QUALITY     Minimum average genotype quality for a variant (default: 10)
  --processes PROCESSES
                        Number of processes. I.e., how many CPUs are available? ~1GB RAM required per process (default: 1)
  --reference_genome {GRCh37,GRCh38}
                        GRCh37/8 (default: None)
  --phased              Use phase information (default: False)
  --microarray          Input is from a microarray. (default: False)
  --debug               Enable debug logging. If not set, logging will be at info level. (default: False)
  --validate            Enable VCF validation. Doubles run time. Might help you identify input issues (default: False)
  --PDFs                Generate a per sample PDF report (default: False)
  --HPAs                Generate results for HPA (default: False)
  --RH                  Generate results for RHD and RHCE. WARNING! EXPERIMENTAL! Based on SNV and small indel
