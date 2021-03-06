﻿Chromothripsis
==============

Simple collection of scripts to analyze data generated by the Cytogenomics Lab
(LIM03/HCFMUSP). At the moment, it covers some bead arrays from Illumina, aCGH
from Agilent and some RNA-Seq data from a Illumina HiSeq. The code is still on
its infancy and I'm not very proficient in R. So, things will look a bit rude
for some time. 

If you find this code useful, don't be shy! Use it! Comments are, of course,
appreciated. 

Licenses will be very liberal, so don't worry about it.

Pipelines
=========

SurePrintG3Custom
--------------------------------

Covers Agilent SurePrint aCGH arrays. It's intended to:
 - Read arrays
 - Preprocessing using limma
 - Genomic wave correction using ArrayTV
 - Segmentation/copy state calling using DNAcopy/CGHcall (not working)
 - Post-processing to get a sheet of regions, state and coverage (partially
 working)

Agilent arrays are particularly susceptible to genomic waves. 

HumanCytoSNP12
--------------

Covers Illumina arrays using crlmm/ArrayTV/VanillaICE. Results are a sheet of
regions, state and coverage:
 - Genotype with arrays crlmm/genotype.Illumina
 - Estimate copy number with crlmmCopynumber
 - Correct for genomic waves with ArrayTV
 - Segment copy state calling with VanillaICE
 - Parse results per state in all samples (coverage in mind)
 - Parse results per sample (dx in mind)

TruSeq
------

Covers RNA-seq made with TruSeq/HiScanSq aimed at variants discovery (SNP,
splicing, etc.) and differential expression:
 - Quality control with ShortRead
 - Mapping with align/Rsubread
 - Splicing junction finding with subjunc/Rsubread 
 - SNP discovery with exactSNP/Rsubread
 - DE with DESeq2

