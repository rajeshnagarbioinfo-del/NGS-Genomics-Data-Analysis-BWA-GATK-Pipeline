#!/bin/bash

# =========================
# WES PIPELINE - 
# =========================

# Step 1: Download SRA file

prefetch SRR28640601 -O data
fasterq-dump SRR28640601 -O data/ -e 6 --split-files

# Step 2: Quality check

fastqc SRR28640601*.fastq -o results/qc/
multiqc results/qc/ -o results/qc/

# Step 3: Trimming

java -jar /usr/share/java/trimmomatic.jar PE -threads 8 \
SRR28640601_1.fastq SRR28640601_2.fastq \
results/trimmed/R1_paired.fq results/trimmed/R1_unpaired.fq \
results/trimmed/R2_paired.fq results/trimmed/R2_unpaired.fq \
ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50

# Step 5: Post-trim QC

fastqc results/trimmed/R1_paired.fq results/trimmed/R2_paired.fq -o results/qc/
multiqc results/qc/ -o results/qc/

# step 6: Reference genome & indexing

mkdir ref
cd ref
wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa human38.fa

#step 7: Indexing

bwa index human38.fa
samtools faidx human38.fa
gatk CreateSequenceDictionary \
-R human38.fa

# Step 4: Alignment
bwa mem -t 8 \
-R "@RG\tID:1\tSM:sample1\tPL:ILLUMINA" \
ref/human38.fa \
results/trimmed/R1_paired.fq R2_paired.fq > results/bam/aligned.sam

# Step 5: BAM processing

samtools view -Sb results/bam/aligned.sam > results/bam/aligned.bam
samtools sort results/bam/aligned.bam -o results/bam/sorted.bam
samtools index results/bam/sorted.bam

# Step 6: Mark duplicates

gatk MarkDuplicates \
-I results/bam/sorted.bam \
-O results/bam/dedup.bam \
-M results/reports/metrics.txt

samtools index results/bam/dedup.bam

# Step 7: Variant Calling
gatk HaplotypeCaller \
-R ref/human38.fa \
-I results/bam/dedup.bam \
-O results/vcf/raw.g.vcf \
-ERC GVCF

# Step 8: Filtering

gatk VariantFiltration \
-R ref/human38.fa \
-V results/vcf/raw.g.vcf \
-O results/vcf/filtered.vcf \
--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
--filter-name "BasicFilter"

# Step 9: Annotation (VEP)
vep -i results/vcf/filtered.vcf \
-o results/annotation/annotated.vcf \
--cache --offline \
--assembly human38
