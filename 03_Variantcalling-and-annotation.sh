#!/usr/bin/bash

# Get user input for necessary paths and parameters
echo -n "Enter Deliverables Directory Path: "
read deliverables_dir

echo -n "Enter BAM data Directory Path: "
read rwDir

# Get file extensions for BAM file
for i in `ls ${rwDir} | grep bqsr_reads.bam`;do echo $i;done
echo -n "Please look into the read names and provide the extension after sample name of BAM Example: dedup_bqsr_reads.bam: "
read "bamextn"


# Get resource allocation details
echo -n "number of threads: "
read cpu

# Create the required directory structure

mkdir -p ${deliverables_dir}/04_Variantcalling
mkdir -p ${deliverables_dir}/00_logs


# Defined paths of Tools & Inputs & Reference

Gatktool="{Path to container tool}/gatk_4.1.3.0.sif"                            # Path to GATK singularity container
FuncotatorDb="{Path to database folder}/funcotator_dataSources.v1.6.20190124s"  # Path to gatk funcotator database
Reference="{Path to Human reference fasta}/HG38/hg38.fa"                        # Path to human reference fasta
logpath="${deliverables_dir}/00_logs"

# Loop through each sample based on the BAM file extension
for i in $(ls $rwDir | grep ${bamextn}); do
    sample=`basename $i ${bamextn}`
    echo "----------------------------------------"
    echo "Processing sample: $sample"
    echo "----------------------------------------"

    bamin=${rwDir}/${sample}${bamextn}
  
    # --- Step 1: Run VariantCalling using HaplotypeCaller ---
    haplotypeok=${logpath}/${sample}_haplotypecaller.ok
    if [ ! -f ${haplotypeok} ]; then
        echo "Running GATK HaplotypeCaller on dedup BAM of ${sample}..."
            ( /usr/bin/time -v ${Gatktool} gatk --java-options "-Xmx8G" HaplotypeCaller -R ${Reference} -I ${bamin} -O ${deliverables_dir}/04_Variantcalling/${sample}_raw_variants.vcf --native-pair-hmm-threads ${cpu} ) 2>&1 | tee ${logpath}/${sample}_haplotypecaller.log && touch ${logpath}/${sample}_haplotypecaller.ok
    else
        echo "Haplotypecaller already completed for ${sample}. Skipping."
    fi
    

    # --- Step 2: Run VariantFiltration on Raw vcf ---
    # This step will only run if the variantcalling step was successful (the .ok file exists)
VariantFiltrationok=${logpath}/${sample}_VariantFiltration.ok
    if [ -f "${haplotypeok}" ] && [ ! -f "${VariantFiltrationok}" ]; then
        echo "Running GATK VariantFiltration on raw vcf ${sample}..."
           ( /usr/bin/time -v ${Gatktool} gatk --java-options "-Xmx8G" VariantFiltration -R ${Reference} -V ${deliverables_dir}/04_Variantcalling/${sample}_raw_variants.vcf -O ${deliverables_dir}/04_Variantcalling/${sample}_filtered_variants.vcf -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 60.0" -filter-name "MQ_filter" -filter "MQ < 40.0" -filter-name "SOR_filter" -filter "SOR > 4.0" -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" -genotype-filter-expression "DP < 10" -genotype-filter-name "DP_filter" -genotype-filter-expression "GQ < 10" -genotype-filter-name "GQ_filter" ) 2>&1 | tee ${logpath}/${sample}_VariantFiltration.log && touch ${logpath}/${sample}_VariantFiltration.ok
    else
        echo "VariantFiltration failed for ${sample}. Skipping VariantFiltration."
    fi
        
        # --- Step 3: VariantCalling using Funcotator on filtered vcf ---
    # This step will only run if the VariantFiltration step was successful (the .ok file exists)
Funcotatorok=${logpath}/${sample}_Funcotator.ok
    if [ -f "${VariantFiltrationok}" ] && [ ! -f "${Funcotatorok}" ]; then
        echo "Running GATK Variant annotation using Funcotator on filter vcf ${sample}..."
           ( /usr/bin/time -v ${Gatktool} gatk --java-options "-Xmx8G" Funcotator --variant ${deliverables_dir}/04_Variantcalling/${sample}_filtered_variants.vcf --reference ${Reference} --ref-version hg38 --data-sources-path ${FuncotatorDb} --output ${deliverables_dir}/04_Variantcalling/${sample}_annotated_variants.vcf --output-file-format VCF ) 2>&1 | tee ${logpath}/${sample}_Funcotator.log && touch ${logpath}/${sample}_Funcotator.ok
    else
        echo "VariantCalling using Funcotator failed for ${sample}. VariantCalling Funcotator."
    fi
        
done

echo "----------------------------------------"
echo "All samples processed."
echo "----------------------------------------"
