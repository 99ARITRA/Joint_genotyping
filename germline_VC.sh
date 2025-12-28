#!/usr/bin/bash
# SOF
# This pipeline is meant to do germline variant calling, filtration and annotation

# ----------------------------------------------- #
## SOURCE MODULE SCRIPTS
# ----------------------------------------------- #
# source ./conda_env.sh

# ----------------------------------------------- #
## TWEAKS 
# ----------------------------------------------- #
BOLD_RED='\033[1;31m'
BOLD_GREEN='\033[1;32m'
BOLD_YELLOW='\033[1;33m'
BOLD_BLUE='\033[1;34m'
BOLD_PURPLE='\033[1;35m'
BOLD_CYAN='\033[1;36m'
NC='\033[0m'

# ----------------------------------------------- #
## INPUT DIR PATHS
# ----------------------------------------------- #
sample_sheet="/mnt/d/Bioinformatics/NGS/GENOMICS/pipelines/bash/germline/germline_samplesheet.csv" # stores sample metadata
fasta="/mnt/d/Bioinformatics/NGS/GENOMICS/refs/genome/chr8.fa" # reference genome
index="/mnt/d/Bioinformatics/NGS/GENOMICS/refs/index/chr8_hg19/chr8" # reference genome index
dict="/mnt/d/Bioinformatics/NGS/GENOMICS/refs/genome/chr8.dict" # sequence dictonary
dbsnp="/mnt/d/Bioinformatics/NGS/GENOMICS/refs/process_files/chr8.dbsnp150.hg19.vcf.gz" # dbSNP annotations VCF file


# ----------------------------------------------- #
## OUTPUR DIR PATHS
# ----------------------------------------------- #
output="/mnt/d/Bioinformatics/NGS/GENOMICS/OP_PROJECT"
bam_data=$output/BAM_DATA
prep_reports=$output/PREPROCESSING_REPORTS
vcf_data=$output/VCF_DATA
logs_dir=$output/LOG_FILES
temp=$output/TEMP # To be deleted at the end of the pipeline

# ----------------------------------------------- #
## CHECK INPUT DIR PATHS AND FILES
# ----------------------------------------------- #
inputFiles=($sample_sheet $fasta $dict $dbsnp)
for file in ${inputFiles[@]}; do
    if [[ ! -f $file ]]; then
        echo -e "${BOLD_RED}>>>[FILE_NOT_FOUND_ERROR] $file is  not present${NC}"
        exit 1
    else
        echo -e "${BOLD_CYAN}>>> $file exist${NC}"
    fi
done
echo -e "\n"
# ----------------------------------------------- #
CHECK_FILE() {
    if [[ ! -f $file ]]; then
        echo -e "${BOLD_RED}>>> $file is  not present${NC}"
        exit 1
    else
        echo -e "${BOLD_CYAN}>>> $file exist${NC}"
    fi
}

# ----------------------------------------------- #
## PIPELINE PARAMETERS 
# ----------------------------------------------- #
# condaLocation=/home/aritra_palodhi/miniforge3/etc/profile.d/conda.sh
# pipelineEnv=ngs
threads=6

# ----------------------------------------------- #
## LOG REPORTS
# ----------------------------------------------- #
LOGS() {
    processLog=$logs_dir/process.log
    $1 2>> $processLog
    echo -e "\n" >> $processLog
}
# ----------------------------------------------- #   
START_TIME() {
    traceLog=$logs_dir/trace.log
    st=$(date +%d-%m-%Y_%H:%M:%S)
    echo "$1 START TIME: $st"
    echo -e "$1\t$st" >>$traceLog
}
# ----------------------------------------------- #
END_TIME() {
    et=$(date +%d-%m-%Y_%H:%M:%S)
    echo "$1 END TIME: $et"
    echo -e "$1\t$et" >>$traceLog
}

# ----------------------------------------------- #
## CREATE THE OUTPUT DIRECTORIES
# ----------------------------------------------- #
BUILD_DIR() {
    mkdir -p $output # Somatic variant calling reports and files output directory
    mkdir -p $bam_data # Storing BAM files
    mkdir -p $prep_reports # Storing preprocessing reports
    mkdir -p $vcf_data # Storing unfiltered and filtered VCF files
    mkdir -p $logs_dir # Log report
    mkdir -p $temp # Directory to store intermediate files in pipeline. Later to be deleted during pipeline run
}

# ----------------------------------------------- #
## STEP-1: MAPPING
# ----------------------------------------------- #
ALIGNMENT() {
    if [[ ! -f $bam_data/${sample}_mapped_sorted.bam || ! -f $bam_data/${sample}_bqsr_sorted.bam ]]; then
        echo -e "${BOLD_BLUE}>>> STEP-1 -->>> Mapping ${NC}\n"
        START_TIME "BWA_MEM"
        bwa mem -t $threads $index $fr $rr  > $bam_data/${sample}.sam
        END_TIME "BWA_MEM"
        # --------------------------------------------------------------- #
        START_TIME "SAMTOOLS"
        samtools view -@ $threads -F 0x4 -h  $bam_data/${sample}.sam | samtools sort -@ $threads > $bam_data/${sample}_mapped_sorted.bam
        # --------------------------------------------------------------- #
        samtools index -@ $threads $bam_data/${sample}_mapped_sorted.bam
        # --------------------------------------------------------------- #
        samtools flagstats -@ $threads $bam_data/${sample}_mapped_sorted.bam > $BAMreports/${sample}_mapped.bam.stats
        END_TIME "SAMTOOLS"
    else
        echo -e "${BOLD_RED}Skipping MAPPING step${NC}\n"
    fi
}

# ----------------------------------------------- #
## STEP-2: BAM PROCESSING
# ----------------------------------------------- #

GATK_ARRG() {
    START_TIME "GATK_ARRG"
    gatk AddOrReplaceReadGroups -I $bam_data/${sample}_mapped_sorted.bam -O $bam_data/${sample}_rg.bam \
        --RGID rg_${sample}   --RGPL illumina  --RGSM ${sample}  --RGPU unit_${sample} --RGLB lib_${sample}
    END_TIME "GATK_ARRG"
}
# ---------------------------------------------------------------------------------------- #
GATK_MD() {
    START_TIME "GATK_MD"
    gatk MarkDuplicates -I $bam_data/${sample}_rg.bam  -M $bam_data/${sample}_dup_metrics.txt  \
        -O $bam_data/${sample}_deduplicated.bam --CREATE_INDEX True --REMOVE_DUPLICATES True
    # ---------------------------------------------------------------------------------------- #
    samtools sort -@ $threads $bam_data/${sample}_deduplicated.bam > $bam_data/${sample}_deduplicated_sorted.bam
    # ---------------------------------------------------------------------------------------- #
    samtools index -@ $threads $bam_data/${sample}_deduplicated_sorted.bam
    # --------------------------------------------------------------- #
    samtools stats -@ $threads $bam_data/${sample}_deduplicated_sorted.bam > $bam_data/${sample}_deduplicated.bam.stats
    END_TIME "GATK_MD"
}
# ------------------------------------------------------------------------------- #
GATK_BQSR() {
    START_TIME "GATK_BQSR"
    gatk BaseRecalibrator -I $bam_data/${sample}_deduplicated_sorted.bam \
            --known-sites $dbsnp \
            -O $BQSRreports/${sample}_BQSR.recalibration.table \
            -R $fasta
     END_TIME "GATK_BQSR"
}
# ---------------------------------------------------------------------------------------- #
GATK_APPLY_BQSR() {
    START_TIME "GATK_APPLY_BQSR"
    gatk ApplyBQSR -R $fasta -I $bam_data/${sample}_deduplicated_sorted.bam \
            --bqsr-recal-file $BQSRreports/${sample}_BQSR.recalibration.table \
            -O $bam_data/${sample}_bqsr.bam
    # ---------------------------------------------------------------------------------------- #
    samtools sort -@ $threads $bam_data/${sample}_bqsr.bam > $bam_data/${sample}_bqsr_sorted.bam
    #  ---------------------------------------------------------------------------------------- #
    samtools index -@ $threads $bam_data/${sample}_bqsr_sorted.bam
    # --------------------------------------------------------------- #
    samtools stats -@ $threads $bam_data/${sample}_bqsr_sorted.bam > $BAMreports/${sample}_bqsr.bam.stats
    END_TIME "GATK_APPLY_BQSR"
    # --------------------------------------------------------------- #
}


BAM_PROCESSING() {
    if [[ ! -f $bam_data/${sample}_deduplicated_sorted.bam  || ! -f $bam_data/${sample}_bqsr_sorted.bam ]]; then
        echo -e "${BOLD_BLUE}>>> STEP-2 -->>> BAM file manipulation ${NC} \n"
        GATK_ARRG # Manipulate BAM records with sample read groups
        GATK_MD # Mark and remove PCR and optical duplicates
        GATK_BQSR # Generate BQSR table
        GATK_APPLY_BQSR # Recalibrate Base quality scores
    else
        echo -e "${BOLD_RED}Skipping BAM PROCESSING step${NC}\n"
    fi
}

# --------------------------------------------------------- #
## MOVE AND REMOVE THE TEMPORARY FILES
# --------------------------------------------------------- #

TEMP_FILES() {
    mv $bam_data/${sample}.sam $temp
    mv $bam_data/${sample}_mapped_sorted.bam $temp
    mv $bam_data/${sample}_mapped_sorted.bam.bai $temp
    mv $bam_data/${sample}_rg.bam $temp
    mv $bam_data/${sample}_deduplicated_sorted.bam $temp
    mv $bam_data/${sample}_deduplicated_sorted.bam.bai $temp
    mv $bam_data/${sample}_bqsr.bam $temp
    mv $vcf_data/${sample}_germline.vcf.gz $temp
    mv $vcf_data/${sample}_germline.vcf.gz.tbi $temp    
    mv $vcf_data/${sample}_germline_normalized.vcf.gz $temp
    mv $vcf_data/${sample}_germline_normalized.vcf.gz.tbi $temp
    mv $vcf_data/${sample}_germline_ft.vcf.gz $temp
    mv $vcf_data/${sample}_germline_ft.vcf.gz.tbi $temp
    # --------------------------------------------------------- #
    rm -r $temp  
}


# --------------------------------------------------------------------------------------------------- #
## RUN ITERATION THROUGH TUMOR AND NORMAL SAMPLES FROM MAPPING TO BAM PROCESSING
# -------------------------------------------------------------------------------------------------- #

PRE_VC() {
    sampleList=$(tail -n +2 $sample_sheet)

    for entry in $sampleList; do
        IFS=',' read -r SAMPLENAME R1 R2 RGID RGPU RGLB PLATFORM <<< $entry
        sample=$SAMPLENAME
        fr=$R1
        rr=$R2
        
        echo -e "${BOLD_BLUE}---------------------------------------------------------------------------------------${NC}\n"
        echo -e "${BOLD_PURPLE} GERMLINE SAMPLE: ${BOLD_GREEN}$sample ${NC}\n"
        echo -e "${BOLD_PURPLE} READ 1: ${BOLD_GREEN}$fr ${NC}\n"
        echo -e "${BOLD_PURPLE} READ 2: ${BOLD_GREEN}$rr ${NC}\n"
        echo -e "${BOLD_BLUE}---------------------------------------------------------------------------------------${NC}\n"

        # ---------------------------------- #
        ## PREPROCESSING PIPELINE
        # ---------------------------------- #
        BUILD_DIR
        # ----------------------------------------------- #
        LOGS "ALIGNMENT"
        # ----------------------------------------------- #
        LOGS "BAM_PROCESSING"
    done
}


# ----------------------------------------------- #
## STEP-3: JOINT GENOTYPING 
# ----------------------------------------------- #
GATK_HAPLOTYPECALLER() {
    if [[ ! -f $vcf_data/${sample}_norm_ft.vcf.gz ]]; then
        START_TIME "GATK_HAPLOTYPECALLER"
        gatk HaplotypeCaller -R $fasta \
                    -I $bam_data/${sample}_bqsr_sorted.bam \
                    -O $vcf_data/${sample}_germline.g.vcf.gz \
                    -ERC GVCF
        END_TIME "GATK_HAPLOTYPECALLER"
    fi
}
# ----------------------------------------------- #
COMBINE_GVCFS() {
    gvcfs=($vcf_data/*_germline.g.vcf)
    # ----------------------------------------------- #
    gatk CombineGVCFs -R $fasta \
                    -V gvcfs[0] \
                    -V gvcfs[1] \
                    -V gvcfs[2] \
                    -O $vcf_data/germline_trio.g.vcf.gz
}
# ----------------------------------------------- #
GENOTYPE_GVCFS() {
    gatk GenotypeGVCFs -R $fasta \
                    -V $vcf_data/germline_trio.g.vcf.gz \
                    -O $vcf_data/joint_genotyped_trio.vcf.gz
    # ----------------------------------------------- #
    gatk IndexFeatureFile -I $vcf_data/joint_genotyped_trio.vcf.gz
}
# ----------------------------------------------- #
NORMALIZE_VCF() {
        # ----------------------------------------------- #
        # STEP-3B: NORMALIZE AND ADD TAGS TO VCF
        # ----------------------------------------------- #
        bcftools norm --threads $threads -c w -f ${fasta} -m-any -o $vcf_data/joint_genotyped_trio_normalized.vcf.gz -Oz $vcf_data/joint_genotyped_trio.vcf.gz
        # ----------------------------------------------- #
        bcftools +fill-tags --threads $threads -o $vcf_data/joint_genotyped_trio_fill-tags.vcf.gz -Oz $vcf_data/joint_genotyped_trio_normalized.vcf.gz
        # ----------------------------------------------- #
        tabix -f --threads $threads -p vcf $vcf_data/joint_genotyped_trio_fill-tags.vcf.gz
        END_TIME "NORMALIZE_VCF"
}    


JOINT_GENOTYPING() {
    echo -e "Step 3 -->>> Germline variant calling: Joint Genotyping"
    GATK_HAPLOTYPECALLER # Call germline variants from each sample
    COMBINE_GVCFS # Combine the sample VCFs
    GENOTYPE_GVCFS # Combined VCF genotyping
    NORMALIZE_VCF # VCF processing

}

# ----------------------------------------------- #
## STEP-4: VARIANT FILTRATION
# ----------------------------------------------- #
# HARD FILTERS FOR HAPLOTYPECALLER OUTPUT
GATK_VARIANT_FILTRATION() {
    START_TIME "GATK_VARIANT_FILTRATION"
    echo -e "${BOLD_YELLOW}>>> STEP 4 -->>> Germline Variant Filtration ${NC}\n"
    gatk VariantFiltration -R $fasta -V $vcf_data/joint_genotyped_trio_fill-tags.vcf.gz -O $vcf_data/joint_genotyped_trio_filtered.vcf.gz \
        --filter-expression  'QUAL < 30.0' --filter-name "Low_Variant_Quality" \
        --filter-expression "MQ < 25.0" --filter-name "Low_Mapping_Quality" \
        --filter-expression "QD < 30.0" --filter-name "Low_Quality_by_Depth" \
        --filter-expression "DP < 1" --filter-name "Low_Read_Depth" \
    END_TIME "GATK_VARIANT_FILTRATION"
}

# ------------------------------------------- #
## CREATE WORKFLOW AND RUN THE PIPELINE
# ------------------------------------------- #

GERMLINE() {
    # NGS_ENV
    echo -e "${BOLD_BLUE}-------------------------------------------------------------------------------------------- "
    echo -e "<<<${BOLD_BLUE}      ${BOLD_PURPLE}* * * ${BOLD_YELLOW}GERMLINE VARIANT ANALYSIS ${BOLD_PURPLE}* * *         ${BOLD_BLUE}>>>"
    echo -e "${BOLD_BLUE}-------------------------------------------------------------------------------------------- ${NC}"                                                     
    PRE_VC # Preprocessing before Variant calling
    LOGS "JOINT_GENOTYPING" # Somatic Variant calling
    LOGS "GATK_VARIANT_FILTRATION" # Somatic variant filtration
    TEMP_FILES
    # ---------------------------------- #
}

GERMLINE

# EOF