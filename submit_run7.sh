#!/bin/sh
 
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=16:mem=32gb
## This tells the batch manager to re-run job with parameter varying from 1 to N in steps on step-size
#PBS -J 1-95
 
## OTHER OPTIONAL PBS DIRECTIVES
 
module load bwa
module load samtools
module load java
module load picard
module load anaconda3/personal
 
/rds/general/project/fisher-aspergillus-analysis/live/bwa_mem_picard_align_run7.sh $(head -$PBS_ARRAY_INDEX /rds/general/project/fisher-aspergillus-analysis/live/arglist_run7.txt | tail -1)
