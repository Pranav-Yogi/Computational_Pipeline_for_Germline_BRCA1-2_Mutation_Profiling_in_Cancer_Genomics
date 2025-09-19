# Computational_Pipeline_for_Germline_BRCA1-2_Mutation_Profiling_in_Cancer_Genomics
This project implements a bioinformatics pipeline for germline BRCA1/2 mutation identification from NGS data. Raw sequencing reads undergo quality assessment with FastQC and preprocessing with Fastp to ensure high-quality input. Reads are then aligned to the human reference genome (hg38) using BWA-MEM, followed by sorting, duplicate marking, and base quality recalibration with Samtools and GATK. Variant discovery is performed with GATK HaplotypeCaller, and raw variants are subjected to stringent filtering to minimize false positives. High-confidence variants are then annotated using GATK Funcotator with reference databases, focusing on clinically relevant BRCA1/2 mutations. The pipeline provides an end-to-end workflow from raw data to annotated germline variants, enabling downstream analysis in cancer risk assessment and precision medicine research.



## EXAMPLE TREE ##
# Above Scripts create folder structure like below 
Deliverables/
├── 00_logs
│   ├── SRR27970273SB_alignment.log
│   ├── SRR27970273SB_alignment.ok
│   ├── SRR27970273SB_Applybqsr.log
│   ├── SRR27970273SB_Applybqsr.ok
│   ├── SRR27970273SB_baserecal.log
│   ├── SRR27970273SB_baserecal.ok
│   ├── SRR27970273SB_buildbamind.log
│   ├── SRR27970273SB_buildbamind.ok
│   ├── SRR27970273SB_fastp.log
│   ├── SRR27970273SB_Funcotator.log
│   ├── SRR27970273SB_Funcotator.ok
│   ├── SRR27970273SB_haplotypecaller.log
│   ├── SRR27970273SB_haplotypecaller.ok
│   ├── SRR27970273SB_MarkDuplicatesSpark.log
│   ├── SRR27970273SB_MarkDuplicatesSpark.ok
│   ├── SRR27970273SB_processed_data_fastqc.log
│   ├── SRR27970273SB_processed_data_fastqc.ok
│   ├── SRR27970273SB_processed_fastp.ok
│   ├── SRR27970273SB_VariantFiltration.log
│   ├── SRR27970273SB_VariantFiltration.ok
│   ├── SRR31058545SB_alignment.log
│   ├── SRR31058545SB_alignment.ok
│   ├── SRR31058545SB_Applybqsr.log
│   ├── SRR31058545SB_Applybqsr.ok
│   ├── SRR31058545SB_baserecal.log
│   ├── SRR31058545SB_baserecal.ok
│   ├── SRR31058545SB_buildbamind.log
│   ├── SRR31058545SB_buildbamind.ok
│   ├── SRR31058545SB_fastp.log
│   ├── SRR31058545SB_Funcotator.log
│   ├── SRR31058545SB_Funcotator.ok
│   ├── SRR31058545SB_haplotypecaller.log
│   ├── SRR31058545SB_haplotypecaller.ok
│   ├── SRR31058545SB_MarkDuplicatesSpark.log
│   ├── SRR31058545SB_MarkDuplicatesSpark.ok
│   ├── SRR31058545SB_processed_data_fastqc.log
│   ├── SRR31058545SB_processed_data_fastqc.ok
│   ├── SRR31058545SB_processed_fastp.ok
│   ├── SRR31058545SB_VariantFiltration.log
│   └── SRR31058545SB_VariantFiltration.ok
├── 00_raw
│   ├── SRR27970273SB_R1.fastq.gz
│   ├── SRR27970273SB_R2.fastq.gz
│   ├── SRR31058545SB_R1.fastq.gz
│   └── SRR31058545SB_R2.fastq.gz
├── 01_QC
│   └── processed_read_qc
│       ├── fastp_reports
│       │   ├── SRR27970273SB_fastp.html
│       │   ├── SRR27970273SB_fastp.json
│       │   ├── SRR31058545SB_fastp.html
│       │   └── SRR31058545SB_fastp.json
│       └── fastqc_reports
│           ├── SRR27970273SB
│           │   ├── SRR27970273SB_R1.fastp_fastqc.html
│           │   ├── SRR27970273SB_R1.fastp_fastqc.zip
│           │   ├── SRR27970273SB_R2.fastp_fastqc.html
│           │   └── SRR27970273SB_R2.fastp_fastqc.zip
│           └── SRR31058545SB
│               ├── SRR31058545SB_R1.fastp_fastqc.html
│               ├── SRR31058545SB_R1.fastp_fastqc.zip
│               ├── SRR31058545SB_R2.fastp_fastqc.html
│               └── SRR31058545SB_R2.fastp_fastqc.zip
├── 02_processed_data
│   ├── SRR27970273SB_R1.fastp.fastq.gz
│   ├── SRR27970273SB_R2.fastp.fastq.gz
│   ├── SRR31058545SB_R1.fastp.fastq.gz
│   └── SRR31058545SB_R2.fastp.fastq.gz
├── 03_Alignment
│   ├── SRR27970273SB_dedup_bqsr_reads.bai
│   ├── SRR27970273SB_dedup_bqsr_reads.bam
│   ├── SRR27970273SB_dedup_reads.bai
│   ├── SRR27970273SB_dedup_reads.bam
│   ├── SRR27970273SB_dedup_reads.bam.bai
│   ├── SRR27970273SB_dedup_reads.bam.sbi
│   ├── SRR27970273SB_recal_data.table
│   ├── SRR27970273SB_sort_reads.bam
│   ├── SRR31058545SB_dedup_bqsr_reads.bai
│   ├── SRR31058545SB_dedup_bqsr_reads.bam
│   ├── SRR31058545SB_dedup_reads.bai
│   ├── SRR31058545SB_dedup_reads.bam
│   ├── SRR31058545SB_dedup_reads.bam.bai
│   ├── SRR31058545SB_dedup_reads.bam.sbi
│   ├── SRR31058545SB_recal_data.table
│   └── SRR31058545SB_sort_reads.bam
└── 04_Variantcalling
    ├── SRR27970273SB_annotated_variants.vcf
    ├── SRR27970273SB_annotated_variants.vcf.idx
    ├── SRR27970273SB_filtered_variants.vcf
    ├── SRR27970273SB_filtered_variants.vcf.idx
    ├── SRR27970273SB_raw_variants.vcf
    ├── SRR27970273SB_raw_variants.vcf.idx
    ├── SRR31058545SB_annotated_variants.vcf
    ├── SRR31058545SB_annotated_variants.vcf.idx
    ├── SRR31058545SB_filtered_variants.vcf
    ├── SRR31058545SB_filtered_variants.vcf.idx
    ├── SRR31058545SB_raw_variants.vcf
    └── SRR31058545SB_raw_variants.vcf.idx

