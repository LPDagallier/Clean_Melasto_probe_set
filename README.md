Cleaning of the Melastomataceae probe set for target enrichment
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

This repository details the cleaning process carried out for the
Melastomataceae probe set.

### Cleaning process

The cleaning process is detailed [here](Probe_set_cleaning_final.md).

### Clean probe set

The new and clean probe set is available here:
[`PROBE_SET_CLEAN.FNA`](CLEAN_PROBE_SET/PROBE_SET_CLEAN.FNA)
(nucleotides version) and
[`PROBE_SET_CLEAN_prot.FAA`](CLEAN_PROBE_SET/PROBE_SET_CLEAN_prot.FAA)
(amino acid version).

**N.B.** The purpose of this clean probe set is to be used
bioinformatically to recover targeted sequences from sequencing reads,
but not to physically target the DNA in vitro.

**Additional note.** It might be interesting to remove short sequences
from the probe set, e.g. with `hybpiper fix_targetfile`
(<https://github.com/mossmatters/HybPiper/wiki/Troubleshooting,-common-issues,-and-recommendations#14-fixing-and-filtering-your-target-file>)

### Comparison between the old and new probe set

See details [here](Comparison_probe_set_old_new.md).

![statplot](Comparison_probe_set_old_new_files/figure-gfm/statplot-1.png)
**Figure 1**. Summary of recovery statistics computed with HybPiper for
the assemblies with the old probe set (blue) and the new probe set in
nucleotide format (yellow), and with the new probe set in amino-acids
format (orange). **A**: number of loci with mapped reads, **B**: number
of loci with assembled sequences, and **C**: number of loci with
assembled sequences equal or longer to 75% of the length of their locus
reference in the probe set. Burrow-Wheeler aligner (bwa) was used to map
the reads with nucleotide probe sets, and Diamond was used for the
amino-acids probe set. Numbers right to the boxplots are the median
value.

### How to cite

> Please cite as: *in prep*
