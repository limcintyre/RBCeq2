RBCeq2 documentation
====================

.. image:: images/Lifeblood-R_Primary_Keyline_RGB.png

.. warning::
   NOT FOR CLINICAL USE

RBCeq2 reads in genomic variant data in the form of variant call files (VCF) and outputs blood group (BG) genotypes and phenotypes. The user docs (Word) explains how RBCeq2 constructs possible allele combinations and then filters them until only possible genotype/phenotype combinations remain. Internal to the tool, all variants are based on human genome release GRCh37 or GRCh38, and in the context of a variant ‘_ref’ means no change to the reference nucleotide for the associated position.

Overview
--------

At the highest level RBCeq2 finds all possible alleles, then filters out those that fail certain logic checks. This allows for an auditable trail of why it has reached a certain result. Every effort has been made to be explicit both in encoding alleles in our database and while writing code. This results in verbose but unambiguous results. Some tools in this space employ Baysian or other likihood methods to filter results – we chose not to do that initially but might add it as an optional feature in the future. Last, some liberties have been taken to standardise syntax and nomenclature across blood groups.

This initial release of RBCeq2 is focused on perfecting the calling of International Society for Blood Transfusion (ISBT) defined BG alleles from simple variants (i.e single nucleotide variant -SNVs/small INDEL) that can be found in standard short read or microarray derived VCFs. Further, it supports the use of long read VCFs (i.e large indels and phased data). However, these features are not as polished.

Bugs
----

This software is extensively tested and accurately reports genotypes/phenotypes based on our inhouse definitions of the ‘correct’ answer, however, there are some examples where the ‘correct’ answer is subjective. These docs are detailed – if you find what you think is a bug in the results from RBCeq2 please take the time to understand if it inline with what we intended or not. We will endeavor to fix any black and white bugs in less than one week. Most of these will be rare variants that are encoded wrong in our database. Further, we value any and all feedback and feature requests.

How To
------

Install via pip (python3.12+) or clone the git repository:

.. code-block:: bash

   pip install RBCeq2

Show help:

.. code-block:: bash

   rbceq2 -h

Usage:

.. code-block:: text

   usage: py main.py --vcf example_multi_sample.vcf.gz --out example --reference_genome GRCh37

   options:
     -h, show this help message and exit
     --vcf       Path to vcf file/s
     --out       Prefix for output files (default: None)
     --allele_depth   Allele depth (default: 10)
     --genotype_quality   Genotype quality (default: 10)
     --processes       Number of processes (default: 1)
                       #More is faster, so long as you have that many CPUs and matched RAM (1:1) available
     --reference_genome GRCh37/8 (default: None)
                       #If your result are all reference/wildtype then you probably got this wrong
     --phased    Use phase information (default: False)
     --microarray   Input is from a microarray. (default: False)
     --debug     Enable debug logging. If not set, logging will be at info level. (default: False)
                       #Turn on to generate output like what is used in the examples in the user documentation
     --validate  Enable VCF validation. Doubles run time (default: False)
                       #Not normally needed
     --PDFs Make single sample PDF reports

Output
------

RBCeq2 generates a ``log.txt``, three programmatically passable TSVs, and a folder with one PDF per sample. The TSVs (one for the genotype and two for phenotype [numeric and alphanumeric separately]) have BGs as column names and sample names (from VCF) as the rows.

1000 Genomes Examples
^^^^^^^^^^^^^^^^^^^^^

Genotypes
"""""""""

+---------+--------------------------------------------------+-------------------------+---------------------+
|         | A4GALT                                           | ABCB6                   | ABCC1               |
+=========+==================================================+=========================+=====================+
| GM18501 | A4GALT*01/A4GALT*02, A4GALT*01/A4GALT*02.02      | ABCB6*01/ABCB6*01W.02   | ABCC1*01/ABCC1*01   |
+---------+--------------------------------------------------+-------------------------+---------------------+
| GM18519 | A4GALT*01.02/A4GALT*01.02                         | ABCB6*01/ABCB6*01W.02   | ABCC1*01/ABCC1*01  |
+---------+--------------------------------------------------+-------------------------+---------------------+
| GM18856 | A4GALT*01/A4GALT*01.02                            | ABCB6*01/ABCB6*01       | ABCC1*01/ABCC1*01  |
+---------+--------------------------------------------------+-------------------------+---------------------+

*Note:* In the genotypes TSV, allele pairs are separated by ``/`` and if multiple pairs are possible these are separated by ``,``.

Phenotype Alphanumeric
""""""""""""""""""""""

+---------+----------------+-------+-------+
|         | A4GALT         | ABCB6 | ABCC1 |
+=========+================+=======+=======+
| GM18501 | P1+,Pk+,(P2+)  | Lan+  | WLF+  |
+---------+----------------+-------+-------+
| GM18519 | P1+,Pk+        | Lan+  | WLF+  |
+---------+----------------+-------+-------+
| GM18856 | P1+,Pk+        | Lan+  | WLF+  |
+---------+----------------+-------+-------+

Phenotype Numeric
"""""""""""""""""

+---------+--------+-------+---------+
|         | A4GALT | ABCB6 | ABCC1   |
+=========+========+=======+=========+
| GM18501 |        |       | ABCC1:1 |
+---------+--------+-------+---------+
| GM18519 |        |       | ABCC1:1 |
+---------+--------+-------+---------+
| GM18856 |        |       | ABCC1:1 |
+---------+--------+-------+---------+

In the phenotype TSVs, antigens are separated by ``,`` so if multiple phenotypes are possible they’re separated by `` | `` (e.g., ``DO:1,2 | DO:1,2,5``). Note, not all BGs have both a numeric and alphanumeric phenotype.

Further details
^^^^^^^^^^^^^^^

Please see the RBCeq2 user documentation Word doc