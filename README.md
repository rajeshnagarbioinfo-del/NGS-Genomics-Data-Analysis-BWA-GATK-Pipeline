# NGS-Genomics-Data-Analysis-BWA-GATK-Pipeline
#  Objective
This project analyzes Whole Exome Sequencing (WES) data of an adolescent female whole blood sample to identify rare and potentially pathogenic variants associated with severe anemia.

#  Dataset information
- SRA ID: SRR28640601  
- Platform: Illumina NovaSeq 6000  
- Type: Paired-end WES  
- Organism: Homo sapiens  
- Reference Genome: GRCh38 (Human Genome Assembly)
  
# Pipeline
FASTQ → Quality Control → Adapter Trimming → Alignment (BWA) → Sorting & Indexing → Duplicate Marking → Base Quality Recalibration → Variant Calling (GATK) → Filtering → Annotation (VEP)

# Tools Used
- linux
- SRA Toolkit
- FastQC / MultiQC  
- Trimmomatic  
- BWA  
- Samtools  
- GATK  
- VEP  

# Key Genes Analyzed
HBB, HBA1, HBA2, TMPRSS6, G6PD

# Output
Filtered and annotated variants linked to anemia-related pathways.

# How to Run this file
bash pipeline.sh
