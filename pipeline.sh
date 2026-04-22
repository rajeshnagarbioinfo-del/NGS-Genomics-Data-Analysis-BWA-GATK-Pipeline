#!/bin/bash


# WES PIPELINE - 

# Step 0: Create folders

mkdir -p data results/qc results/multiqc results/trimmed results/bam results/vcf results/annotation results/reports ref

# Step 1: Download SRA file

prefetch SRR28640601 -O data
fasterq-dump SRR28640601 -O data/ -e 6 --split-files

cd data

# Step 2: Quality check

fastqc SRR28640601*.fastq -o results/qc/
multiqc results/qc/ -o results/multiqc

# Step 3: Trimming

java -jar /usr/share/java/trimmomatic.jar PE -threads 8 \
SRR28640601_1.fastq SRR28640601_2.fastq \
results/trimmed/R1_paired.fq results/trimmed/R1_unpaired.fq \
results/trimmed/R2_paired.fq results/trimmed/R2_unpaired.fq \
ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50

# Step 4: Post-trim QC

fastqc results/trimmed/R1_paired.fq results/trimmed/R2_paired.fq -o results/qc/
multiqc results/qc/ -o results/multiqc

cd ..

# step 5: Reference genome & indexing

cd ref
wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa human38.fa

#step 6: Indexing

bwa index human38.fa
samtools faidx human38.fa
gatk CreateSequenceDictionary \
-R human38.fa

cd ..

# Step 7: Alignment

bwa mem -t 8 \
-R "@RG\tID:1\tSM:sample1\tPL:ILLUMINA" \
ref/human38.fa \
results/trimmed/R1_paired.fq results/trimmed/R2_paired.fq > results/bam/aligned.sam

# Step 8: BAM processing

samtools view -Sb results/bam/aligned.sam > results/bam/aligned.bam
samtools sort results/bam/aligned.bam -o results/bam/sorted.bam
samtools index results/bam/sorted.bam

# Step 9: Mark duplicates

gatk MarkDuplicates \
-I results/bam/sorted.bam \
-O results/bam/dedup.bam \
-M results/reports/metrics.txt

samtools index results/bam/dedup.bam

# Step 10: Variant Calling
gatk HaplotypeCaller \
-R ref/human38.fa \
-I results/bam/dedup.bam \
-O results/vcf/raw.g.vcf \
-ERC GVCF

 gatk GenotypeGVCFs \
-R ref/human38.fa \
-V results/vcf/raw.g.vcf \
-O results/vcf/raw.vcf

# Step 11: Filtering

gatk VariantFiltration \
-R ref/human38.fa \
-V results/vcf/raw.vcf \
-O results/vcf/filtered.vcf \
--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
--filter-name "BasicFilter"

bgzip -c results/vcf/filtered.vcf > results/vcf/filtered.vcf.gz
tabix -p vcf results/vcf/filtered.vcf.gz

# step 12 : SNP & INDEL split
#SNP
gatk SelectVariants \
-V filtered.vcf.gz \
-select-type SNP \
-O snp.vcf.gz

tabix -p vcf results/vcf/snp.vcf.gz

#INDEL 
gatk SelectVariants \
-V filtered.vcf.gz \
-select-type INDEL \
-O indel.vcf.gz

tabix -p vcf results/vcf/indel.vcf.gz


# Step 13: SNP Filtering
---------------------------------------------
gatk VariantFiltration \
-V results/vcf/snp.vcf.gz \
-O results/vcf/snp_filtered.vcf.gz \
--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0" \
--filter-name "SNP_Filter"

gatk SelectVariants \
-V results/vcf/snp_filtered.vcf.gz \
--exclude-filtered \
-O results/vcf/snp_pass.vcf.gz


# Step 14: INDEL Filtering
---------------------------------------------
gatk VariantFiltration \
-V results/vcf/indel.vcf.gz \
-O results/vcf/indel_filtered.vcf.gz \
--filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
--filter-name "INDEL_Filter"

gatk SelectVariants \
-V results/vcf/indel_filtered.vcf.gz \
--exclude-filtered \
-O results/vcf/indel_pass.vcf.gz


# Step 15: SNP Annotation (VEP)
-----------------------------------------------
cd ~/dna_pipeline/ensembl-vep

./vep -i ~/dna_pipeline/results/vcf/snp_pass.vcf.gz \
-o ~/dna_pipeline/results/annotation/snp_annotated.vcf \
--cache \
--offline \
--species homo_sapiens \
--assembly GRCh38 \
--vcf \
--symbol \
--hgvs \
--protein \
--sift b \
--polyphen b \
--canonical \
--biotype \
--af \
--fasta ~/dna_pipeline/ref/human38.fa \
--fork 4 \
--force_overwrite


# Step 16: INDEL Annotation (VEP)

./vep -i ~/dna_pipeline/results/vcf/indel_pass.vcf.gz \
-o ~/dna_pipeline/results/annotation/indel_annotated.vcf \
--cache \
--offline \
--species homo_sapiens \
--assembly GRCh38 \
--vcf \
--symbol \
--hgvs \
--protein \
--sift b \
--polyphen b \
--canonical \
--biotype \
--af \
--fasta ~/dna_pipeline/ref/human38.fa \
--fork 4 \
--force_overwrite
