#!/usr/bin/bash

# Get user input for necessary paths and parameters
echo -n "Enter Deliverables Directory Path: "
read deliverables_dir

echo -n "Enter Processed data Directory Path: "
read rwDir

# Get file extensions for R1 and R2 reads
for i in `ls ${rwDir} | grep _R1`;do echo $i;done
echo -n "Please look into the read names and provide the extension after sample name Example: _R1_001.fastq.gz: "
read "R1_extn"

for i in `ls ${rwDir} | grep _R2`;do echo $i;done
echo -n "Please look into the read names and provide the extension after sample name Example: _R2_001.fastq.gz: "
read "R2_extn"

# Get resource allocation details
echo -n "number of threads: "
read cpu

# Defined paths of Tools & Inputs & Reference

bwatool="{Provide your path of container}/bwa_0.7.18--he4a0461_1.sif"                          # Path of singularity container
SamTools="{Provide your path of container}/samtools_1.21--h50ea8bc_0.sif"                      # Path of singularity container
Gatktool="{Provide your path of container}/gatk_4.1.3.0.sif"                                   # Path of singularity container
Reference="{Path of reference}/hg38.fa"                                                        # Path of Indexed Hg38 Human reference genome fasta
Knownsites_dbsnp="{Path of reference}/dbsnp_ucsc_general/Homo_sapiens_assembly38.dbsnp138.vcf" # Path of Known sites Ex. dbSNP UCSC known site vcf
Knownsites_mills1000G="{Path of reference}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"   # Path of Another Knownsite vcf Ex. Mills & 1000G goldstandard indels 

# Create the required directory structure

mkdir -p ${deliverables_dir}/03_Alignment
mkdir -p ${deliverables_dir}/00_logs

# Define log path and adapter path
logpath="${deliverables_dir}/00_logs"

# Loop through each sample based on the R1 file extension
for i in $(ls $rwDir | grep ${R1_extn}); do
    sample=`basename $i ${R1_extn}`
    echo "----------------------------------------"
    echo "Processing sample: $sample"
    echo "----------------------------------------"

    # Define input and output file paths for the current sample
    R1_in=${rwDir}/${sample}${R1_extn}
    R2_in=${rwDir}/${sample}${R2_extn}


    # --- Step 1: Run BWA MEM for alignment AND samtools sort ---
    
    alignmentok=${logpath}/${sample}_alignment.ok
    if [ ! -f ${alignmentok} ]; then
        echo "Running BWA MEM/ Alignment ${sample}..."
            ( /usr/bin/time -v ${bwatool} bwa mem -t ${cpu} -M -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA\tLB:${sample}\tPU:${sample}" ${Reference} ${R1_in} ${R2_in} | ${SamTools} samtools sort -@ 2 -o $deliverables_dir/03_Alignment/${sample}_sort_reads.bam ) 2>&1 | tee ${logpath}/${sample}_alignment.log && touch ${logpath}/${sample}_alignment.ok
    else
        echo "alignment already completed for ${sample}. Skipping."
    fi
    

    # --- Step 2: Run Markduplicatespark on dedup bam ---
    # This step will only run if the alignment step was successful (the .ok file exists)
markdup_ok=${logpath}/${sample}_MarkDuplicatesSpark.ok
    if [ -f "${alignmentok}" ] && [ ! -f "${markdup_ok}" ]; then
        echo "Running GATK markduplicatespark on dedup bam ${sample}..."
           ( /usr/bin/time -v ${Gatktool} gatk MarkDuplicatesSpark -I $deliverables_dir/03_Alignment/${sample}_sort_reads.bam -O $deliverables_dir/03_Alignment/${sample}_dedup_reads.bam ) 2>&1 | tee ${logpath}/${sample}_MarkDuplicatesSpark.log && touch ${logpath}/${sample}_MarkDuplicatesSpark.ok
    else
        echo "Markduplicatespark failed for ${sample}. Skipping markduplicatespark."
    fi
    
    
        # --- Step 3: Run Bam Index on dedup bam ---
    # This step will only run if the Markduplicates spark step was successful (the .ok file exists)
buildbamind_ok=${logpath}/${sample}_buildbamind.ok
    if [ -f "${markdup_ok}" ] && [ ! -f "${buildbamind_ok}" ]; then
        echo "Running GATK markduplicatespark on dedup bam ${sample}..."
           ( /usr/bin/time -v ${Gatktool} gatk BuildBamIndex -I $deliverables_dir/03_Alignment/${sample}_dedup_reads.bam ) 2>&1 | tee ${logpath}/${sample}_buildbamind.log && touch ${logpath}/${sample}_buildbamind.ok
    else
        echo "Build Bam index failed for ${sample}. Build Bam index."
    fi
    

        # --- Step 4: Run BaseRecalibrator on dedup bam ---
    # This step will only run if the buildbamindex step was successful (the .ok file exists)
baserecal_ok=${logpath}/${sample}_baserecal.ok
    if [ -f "${buildbamind_ok}" ] && [ ! -f "${baserecal_ok}" ]; then
        echo "Running GATK markduplicatespark on dedup bam ${sample}..."
           ( /usr/bin/time -v ${Gatktool} gatk --java-options "-Xmx10G" BaseRecalibrator -I $deliverables_dir/03_Alignment/${sample}_dedup_reads.bam -R ${Reference} --known-sites ${Knownsites_dbsnp} --known-sites ${Knownsites_mills1000G} -O $deliverables_dir/03_Alignment/${sample}_recal_data.table ) 2>&1 | tee ${logpath}/${sample}_baserecal.log && touch ${logpath}/${sample}_baserecal.ok
    else
        echo "BaseRecalibrator failed for ${sample}. BaseRecalibrator."
    fi
    

        # --- Step 5: Run ApplyBQSR on dedup bam ---
    # This step will only run if the BaseRecalibrator step was successful (the .ok file exists)
Applybqsr_ok=${logpath}/${sample}_Applybqsr.ok
    if [ -f "${baserecal_ok}" ] && [ ! -f "${Applybqsr_ok}" ]; then
        echo "Running GATK markduplicatespark on dedup bam ${sample}..."
           ( /usr/bin/time -v ${Gatktool} gatk ApplyBQSR -I $deliverables_dir/03_Alignment/${sample}_dedup_reads.bam -R ${Reference} --bqsr-recal-file $deliverables_dir/03_Alignment/${sample}_recal_data.table -O $deliverables_dir/03_Alignment/${sample}_dedup_bqsr_reads.bam ) 2>&1 | tee ${logpath}/${sample}_Applybqsr.log && touch ${logpath}/${sample}_Applybqsr.ok
    else
        echo "ApplyBQSR failed for ${sample}. Applybqsr."
    fi
    
done

echo "----------------------------------------"
echo "All samples processed."
echo "----------------------------------------"
