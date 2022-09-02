#!/bin/sh

PATH=$PATH:/apps/bwa/0.7.15/bin/bwa:/apps/samtools/1.3.1/bin/samtools:/apps/picard/2.6.0/picard.jar:/apps/java/jdk-8u66/bin/java

## SCRIPT TO ALIGN PAIRED END ILLUMINA READS TO A REFERENCE ASSEMBLY USING BQSR WITH KNOWN SITES
## INPUT FILES REQUIRED
## $1 --> Left/1st of paired read in fastq format
## $2 --> Right/2nd of paired read in fastq format
## $3 --> Output prefix
## Versions used:
## BWA aligner: 0.7.15
## SAMTools: 1.2.1
## Picard: 2.6.0
## GATK: 4.2.6.1
#####################################################################

mkdir -p -v /rds/general/project/fisher-aspergillus-results/ephemeral/$3

reference_dir=/rds/general/project/fisher-aspergillus-reference/live
reference=$reference_dir/GCF_000002655.1_ASM265v1_genomic.fa
output=/rds/general/project/fisher-aspergillus-results/ephemeral/$3
results_dir=/rds/general/project/fisher-aspergillus-results/live/Run7/$3
tmp=/rds/general/project/fisher-aspergillus-results/ephemeral/
rgid=A00478-$3-168
rgpu=168

## copy files over to local scratch

cp -v $1 $output/R1.fastq.gz
cp -v $2 $output/R2.fastq.gz

bwa mem -M $reference $output/R1.fastq.gz $output/R2.fastq.gz > $output/$3.sam

samtools import $reference_dir/GCF_000002655.1_ASM265v1_genomic.fa.fai $output/$3.sam $output/$3.bam

samtools sort $output/$3.bam $output/$3.sorted
samtools index $output/$3.sorted.bam

picard AddOrReplaceReadGroups INPUT=$output/$3.sorted.bam OUTPUT=$output/$3.fixed.bam SORT_ORDER=coordinate RGID=$rgid RGLB=dnaseq RGPL=illumina RGSM=WGS RGPU=$rgpu CREATE_INDEX=TRUE VALIDATION_STRINGENCY=SILENT

picard MarkDuplicates INPUT=$output/$3.fixed.bam OUTPUT=$output/$3.sorted.marked.bam METRICS_FILE=$output/picard_info.txt REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=TRUE

#Remove extra files
rm $output/$3.sam
rm $output/$3.bam 

# Activate GATK
source activate gatk4_env

#Call variants
gatk HaplotypeCaller --tmp-dir $tmp -R $reference -I $output/$3.sorted.marked.bam -O $output/$3.raw_variants.vcf --pcr-indel-model NONE -ploidy 1 -stand-call-conf 30 -mbq 20 -A QualByDepth -XL $reference_dir/GCF_000002655.1_ASM265v1_genomic.repeat.intervals 

#Extract SNPs and INDELs
gatk SelectVariants --tmp-dir $tmp -R $reference -V $output/$3.raw_variants.vcf --select-type-to-include SNP -O $output/$3.raw_snps.vcf 

gatk SelectVariants --tmp-dir $tmp -R $reference -V $output/$3.raw_variants.vcf --select-type-to-include INDEL -O $output/$3.raw_indels.vcf 

#Filter SNPs and INDELs
gatk VariantFiltration --tmp-dir $tmp -R $reference -V $output/$3.raw_snps.vcf --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0" --filter-name LowConf -O $output/$3.filtered_snps.vcf 

gatk VariantFiltration --tmp-dir $tmp -R $reference -V $output/$3.raw_indels.vcf --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filter-name LowConf -O $output/$3.filtered_indels.vcf 

#Base Quality Score Recalibration (BQSR) #1
gatk BaseRecalibrator --tmp-dir $tmp -R $reference -I $output/$3.sorted.marked.bam --known-sites $output/$3.filtered_snps.vcf --known-sites $output/$3.filtered_indels.vcf -O $output/$3.recal_data.table 

#Apply BQSR #1
gatk ApplyBQSR --tmp-dir $tmp -R $reference -I $output/$3.sorted.marked.bam --bqsr-recal-file $output/$3.recal_data.table -O $output/$3.recal_reads.bam 

#Base Quality Score Recalibration (BQSR) #2
gatk BaseRecalibrator --tmp-dir $tmp -R $reference -I $output/$3.recal_reads.bam --known-sites $output/$3.filtered_snps.vcf --known-sites $output/$3.filtered_indels.vcf -O $output/$3.post_recal_data.table 

#Apply BQSR
gatk ApplyBQSR --tmp-dir $tmp -R $reference -I $output/$3.recal_reads.bam -O $output/$3.post_recal_reads.bam --bqsr-recal-file $output/$3.post_recal_data.table 

#Remove extra files
rm $output/$3.sorted.marked.ba*

#Call variants (again)
gatk HaplotypeCaller --tmp-dir $tmp -R $reference -I $output/$3.post_recal_reads.bam -O $output/$3.raw_variants_recal.vcf -ERC GVCF --pcr-indel-model NONE -ploidy 1 -stand-call-conf 30 -mbq 20 -A QualByDepth -XL $reference_dir/GCF_000002655.1_ASM265v1_genomic.repeat.intervals  

#Genotype gVCF
gatk GenotypeGVCFs --tmp-dir $tmp -R $reference -V $output/$3.raw_variants_recal.vcf -O $output/$3.genotyped_variants_recal.vcf

#Extract SNPs for stats purposes
gatk SelectVariants --tmp-dir $tmp -R $reference -V $output/$3.genotyped_variants_recal.vcf --select-type-to-include SNP -O $output/$3.raw_snps_stats.vcf

#Extract SNPs and INDELs
gatk SelectVariants --tmp-dir $tmp -R $reference -V $output/$3.genotyped_variants_recal.vcf --select-type-to-include SNP -O $output/$3.raw_snps_recal.vcf -select 'vc.getGenotype("WGS").getAD().1*1.0 / vc.getGenotype("WGS").getDP() > 0.90'

gatk SelectVariants --tmp-dir $tmp -R $reference -V $output/$3.genotyped_variants_recal.vcf --select-type-to-include INDEL -O $output/$3.raw_indels_recal.vcf 

#Filter SNPs
gatk VariantFiltration --tmp-dir $tmp -R $reference -V $output/$3.raw_snps_recal.vcf -filter "QD < 2.0" --filter-name "LowConf" -filter "FS > 60.0" --filter-name "LowConf" -filter "MQ < 40.0" --filter-name "LowConf" -filter "MQRankSum < -12.5" --filter-name "LowConf" -filter "ReadPosRankSum < -8.0" --filter-name "LowConf" -filter "SOR > 4.0" --filter-name "LowConf" -filter "DP < 5" --filter-name "LowConf" -G-filter "GQ < 50" -G-filter-name "FILTER_GQ-50" -O $output/$3.filtered_snps_final.vcf 

#Filter INDELs
gatk VariantFiltration --tmp-dir $tmp -R $reference -V $output/$3.raw_indels_recal.vcf -filter "QD < 2.0" --filter-name "LowConf" -filter "FS > 200.0" --filter-name "LowConf" -filter "ReadPosRankSum < -20.0" --filter-name "LowConf" -filter "SOR > 10.0" --filter-name "LowConf" -O $output/$3.filtered_indels_final.vcf 

grep PASS $output/$3.filtered_snps_final.vcf | awk '$4=="A"||$4=="C"||$4=="G"||$4=="T"' | awk '$5=="A"||$5=="C"||$5=="G"||$5=="T"' > $output/$3.final_snps.body 
grep "#" $output/$3.filtered_snps_final.vcf > $output/$3.final.head
cat $output/$3.final.head $output/$3.final_snps.body > $output/$3.final_snps.vcf

grep PASS $output/$3.filtered_indels_final.vcf | awk '$4=="A"||$4=="C"||$4=="G"||$4=="T"' | awk '$5=="A"||$5=="C"||$5=="G"||$5=="T"' > $output/$3.final_indels.body
grep "#" $output/$3.filtered_indels_final.vcf > $output/$3.final_indels.head
cat $output/$3.final_indels.head $output/$3.final_indels.body > $output/$3.final_indels.vcf

#Calculate depth of coverage
picard CollectWgsMetrics I=$output/$3.post_recal_reads.bam O=$results_dir/$3.metrics.txt R=$reference

#Collect mapping statistics
samtools flagstat $output/$3.post_recal_reads.bam > $output/flagstat

cp -v $output/$3.genotyped_variants_recal.vcf $results_dir
cp -v $output/$3.post_recal_reads.bam $results_dir
cp -v $output/$3.post_recal_reads.bai $results_dir
cp -v $output/$3.raw_variants_recal.vcf $results_dir
cp -v $output/$3.raw_snps_recal.vcf $results_dir
cp -v $output/$3.raw_indels_recal.vcf $results_dir
cp -v $output/$3.final_indels.vcf $results_dir
cp -v $output/$3.final_snps.vcf $results_dir
cp -v $output/flagstat $results_dir
cp -v $output/$3.filtered_snps_final.vcf $results_dir
cp -v $output/$3.raw_snps_stats.vcf $results_dir

exit 0
