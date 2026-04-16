#!/bin/bash

# =========================
# WES PIPELINE - GRCh38
# =========================

# Step 1: Download
prefetch SRR28640601 -O data/raw/
fasterq-dump SRR28640601 -O data/raw/ -e 6 --split-files
gzip data/raw/*.fastq

# Step 2: QC
fastqc data/raw/*.fastq.gz -o results/qc/
multiqc results/qc/ -o results/qc/

# Step 3: Trimming
trimmomatic PE -threads 6 \
data/raw/*_1.fastq.gz data/raw/*_2.fastq.gz \
data/trimmed/R1_paired.fq.gz data/trimmed/R1_unpaired.fq.gz \
data/trimmed/R2_paired.fq.gz data/trimmed/R2_unpaired.fq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50

# Step 4: Alignment
bwa mem -t 8 reference/GRCh38.fa \
data/trimmed/R1_paired.fq.gz data/trimmed/R2_paired.fq.gz > results/bam/aligned.sam

# Step 5: BAM processing
samtools view -Sb results/bam/aligned.sam > results/bam/aligned.bam
samtools sort results/bam/aligned.bam -o results/bam/sorted.bam
samtools index results/bam/sorted.bam

# Step 6: Mark duplicates
gatk MarkDuplicates \
-I results/bam/sorted.bam \
-O results/bam/dedup.bam \
-M reports/metrics.txt

samtools index results/bam/dedup.bam

# Step 7: Variant Calling
gatk HaplotypeCaller \
-R reference/GRCh38.fa \
-I results/bam/dedup.bam \
-O results/vcf/raw.vcf \
-ERC GVCF

# Step 8: Filtering
gatk VariantFiltration \
-R reference/GRCh38.fa \
-V results/vcf/raw.vcf \
-O results/vcf/filtered.vcf \
--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
--filter-name "BasicFilter"

# Step 9: Annotation (VEP)
vep -i results/vcf/filtered.vcf \
-o results/annotation/annotated.vcf \
--cache --offline \
--assembly GRCh38
