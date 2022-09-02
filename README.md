# WGSpipeline
Pipeline to align paired end short reads to reference genome. Assumes reference genome has already been appropriately indexed.

Inputs:

1 - left/1st of paired read in fastq format
2 - right/2nd of paired read in fastq format
3 - output prefix

You can acquire the RGID and RGPU from the fastq header. By including this now in the alignment, you avoid the need to redo read groups later on (this is important for population analyses such as Fst, tajima's D, nucleotide diversity etc).

Creates final BAM file, calls variants (SNPs and INDELs), coverage file and mapping statistics

This is the basic pipeline, but this is an example of a batch submission. Use with 'arglist_run7.txt' and 'submit_run7.sh' to run a batch submission of many jobs.
