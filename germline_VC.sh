#!/usr/bin/bash

# SOF
# This pipeline is meant to do germline variant calling, filtration and annotation

# ----------------------------------------------- #
## SOURCE MODULE SCRIPTS
# ----------------------------------------------- #

source ./conda_env.sh


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

sampleFile=$PWD/germline_samplesheet.csv # stores sample sample, sample type, paired end file paths
refFasta=$PWD/../../../../ref_files/chr1120.fa # reference genome
refIndex=$PWD/../../../../index/chr1120/chr1120  # reference genome index
seqDict=$PWD/../../../../ref_files/chr1120.dict # sequence dictonary
dbsnp=$PWD/../../../../../../Files/vcf/chr1120_dbsnp138.vcf.gz # dbSNP annotations VCF file


# ----------------------------------------------- #
## OUTPUR DIR PATHS
# ----------------------------------------------- #

output=$PWD/../GERMLINE_OUTPUT
BAMdir=$output/BAMS
BAMreports=$output/BAM_REPORT
BQSRreports=$output/BQSR_REPORTS
VCFdir=$output/VCFS
logDir=$output/LOGS
temp=$output/TEMP # To be deleted at the end of the pipeline


# ----------------------------------------------- #
## CHECK INPUT DIR PATHS AND FILES
# ----------------------------------------------- #

inputFiles=($sampleFile $refFasta $seqDict $dbsnp)
for file in ${inputFiles[@]}; do
    if [[ ! -f $file ]]; then
        echo -e "${BOLD_RED}>>>[FILE_NOT_FOUND_ERROR] $(realpath $file) is  not present${NC}"
        exit 1
    else
        echo -e "${BOLD_CYAN}>>> $(realpath $file) exist${NC}"
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

condaLocation=/home/aritra_palodhi/miniforge3/etc/profile.d/conda.sh
pipelineEnv=ngs
threads=6


# ----------------------------------------------- #
## LOG REPORTS
# ----------------------------------------------- #

LOGS() {
    processLog=$logDir/process.log
    $1 2>> $processLog
    echo -e "\n" >> $processLog
}
# ----------------------------------------------- #   
START_TIME() {
    traceLog=$logDir/trace.log
    st=$(date +%Y-%m-%d_%H:%M:%S)
    echo "$1 START TIME: $st"
    echo -e "$1\t$st" >>$traceLog
}
# ----------------------------------------------- #
END_TIME() {
    et=$(date +%Y-%m-%d_%H:%M:%S)
    echo "$1 END TIME: $et"
    echo -e "$1\t$et" >>$traceLog
}

# ----------------------------------------------- #
## CREATE THE OUTPUT DIRECTORIES
# ----------------------------------------------- #

BUILD_DIR() {
    mkdir -p $output # Somatic variant calling reports and files output directory
    mkdir -p $BAMdir # Storing BAM files
    mkdir -p $BAMreports
    mkdir -p $VCFdir # Storing unfiltered and filtered VCF files
    mkdir -p $BQSRreports # Storing BQSR reports
    mkdir -p $logDir # Log report
    mkdir -p $temp # Directory to store intermediate files in pipeline. Later to be deleted during pipeline run
}


# ----------------------------------------------- #
## STEP-1: MAPPING
# ----------------------------------------------- #

ALIGNMENT() {
    if [[ ! -f $BAMdir/${sample}_mapped_sorted.bam || ! -f $BAMdir/${sample}_bqsr_sorted.bam ]]; then
        echo -e "${BOLD_BLUE}>>> STEP-1 -->>> Mapping ${NC}\n"
        START_TIME "BWA_MEM"
        bwa mem -t $threads $refIndex $fr $rr  > $BAMdir/${sample}.sam
        END_TIME "BWA_MEM"
        # --------------------------------------------------------------- #
        START_TIME "SAMTOOLS"
        samtools view -@ $threads -F 0x4 -h  $BAMdir/${sample}.sam | samtools sort -@ $threads > $BAMdir/${sample}_mapped_sorted.bam
        # --------------------------------------------------------------- #
        samtools index -@ $threads $BAMdir/${sample}_mapped_sorted.bam
        # --------------------------------------------------------------- #
        samtools flagstats -@ $threads $BAMdir/${sample}_mapped_sorted.bam > $BAMreports/${sample}_mapped.flagstats
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
    gatk AddOrReplaceReadGroups -I $BAMdir/${sample}_mapped_sorted.bam -O $BAMdir/${sample}_rg.bam \
        --RGID 10 --RGPL ILLUMINA  --RGSM ${sample}  --RGPU unit10 --RGLB lib10
    END_TIME "GATK_ARRG"
}
# ---------------------------------------------------------------------------------------- #
GATK_MD() {
    START_TIME "GATK_MD"
    gatk MarkDuplicates -I $BAMdir/${sample}_rg.bam  -M $BAMdir/${sample}_dup_metrics.txt  \
        -O $BAMdir/${sample}_rg_markdup.bam --CREATE_INDEX True --REMOVE_DUPLICATES True
    # ---------------------------------------------------------------------------------------- #
    samtools sort -@ $threads $BAMdir/${sample}_rg_markdup.bam > $BAMdir/${sample}_sorted_rg_markdup.bam
    # ---------------------------------------------------------------------------------------- #
    samtools index -@ $threads $BAMdir/${sample}_sorted_rg_markdup.bam
    # --------------------------------------------------------------- #
    samtools flagstats -@ $threads $BAMdir/${sample}_sorted_rg_markdup.bam > $BAMdir/${sample}_markdup.flagstats
    END_TIME "GATK_MD"
}
# ------------------------------------------------------------------------------- #
GATK_BQSR() {
    START_TIME "GATK_BQSR"
    gatk BaseRecalibrator -I $BAMdir/${sample}_sorted_rg_markdup.bam \
            --known-sites $dbsnp \
            -O $BQSRreports/${sample}_BQSR_recalibration.table \
            -R $refFasta
     END_TIME "GATK_BQSR"
}
# ---------------------------------------------------------------------------------------- #
GATK_APPLY_BQSR() {
    START_TIME "GATK_APPLY_BQSR"
    gatk ApplyBQSR -R $refFasta -I $BAMdir/${sample}_sorted_rg_markdup.bam \
            --bqsr-recal-file $BQSRreports/${sample}_BQSR_recalibration.table \
            -O $BAMdir/${sample}_bqsr.bam
    # ---------------------------------------------------------------------------------------- #
    samtools sort -@ $threads $BAMdir/${sample}_bqsr.bam > $BAMdir/${sample}_sorted_bqsr.bam
    #  ---------------------------------------------------------------------------------------- #
    samtools index -@ $threads $BAMdir/${sample}_sorted_bqsr.bam
    # --------------------------------------------------------------- #
    samtools flagstats -@ $threads $BAMdir/${sample}_sorted_bqsr.bam > $BAMreports/${sample}_sorted_bqsr.flagstats
    END_TIME "GATK_APPLY_BQSR"
    # --------------------------------------------------------------- #
}


BAM_PROCESSING() {
    if [[ ! -f $BAMdir/${sample}_sorted_rg_markdup.bam  || ! -f $BAMdir/${sample}_sorted_rg_markdup_bqsr.bam ]]; then
        echo -e "${BOLD_BLUE}>>> STEP-2 -->>> BAM file manipulation ${NC} \n"
        GATK_ARRG # GATK AddorReplaceReadGroups
        GATK_MD # GATK Markduplicates
        GATK_BQSR # GATK BaseQualityScoreRecalibration
        GATK_APPLY_BQSR # GATK ApplyBQSR
    else
        echo -e "${BOLD_RED}Skipping BAM PROCESSING step${NC}\n"
    fi
}

# --------------------------------------------------------- #
## MOVE AND REMOVE THE TEMPORARY FILES
# --------------------------------------------------------- #

TEMP_FILES() {
    mv $BAMdir/${sample}.sam $temp
    mv $BAMdir/${sample}_mapped_sorted.bam $temp
    mv $BAMdir/${sample}_mapped_sorted.bam.bai $temp
    mv $BAMdir/${sample}_mapped_sorted_rg.bam $temp
    mv $BAMdir/${sample}_sorted_rg.bam $temp
    mv $BAMdir/${sample}_sorted_rg_markdup.bam $temp
    mv $BAMdir/${sample}_sorted_rg_markdup.bam.bai $temp
    mv $BAMdir/${sample}_sorted_rg_markdup_bqsr.bam $temp
    mv $VCFdir/${sample}.vcf.gz $temp
    mv $VCFdir/${sample}.vcf.gz.tbi $temp    
    mv $VCFdir/${sample}_norm.vcf.gz $temp
    # --------------------------------------------------------- #
    rm -r $temp  
}


# --------------------------------------------------------------------------------------------------- #
## RUN ITERATION THROUGH TUMOR AND NORMAL SAMPLES FROM MAPPING TO BAM PROCESSING
# -------------------------------------------------------------------------------------------------- #

PRE_VC() {
    sampleList=$(tail -n +2 $sampleFile)

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
## STEP-3: VARIANT CALLING
# ----------------------------------------------- #
# ----------------------------------------------- #
## STEP-3A: VARIANT CALLING USING HAPLOTYPECALLER
# ----------------------------------------------- #
GATK_HAPLOTYPECALLER() {
    if [[ ! -f $VCFdir/${sample}_norm_ft_germline.vcf.gz ]]; then
        echo -e "${BOLD_YELLOW}>>> STEP 3 -->>> Germline Variant Calling ${NC}\n"
        START_TIME "GATK_HAPLOTYPECALLER"
        read -p "Select the type of germline variant calling [Single / Joint genotyping]: " choice
        gatk HaplotypeCaller -R $refFasta \
                    -I $BAMdir/${sample}_sorted_bqsr.bam \
                    -O $VCFdir/${sample}_germline.vcf
        # ----------------------------------------------- #
        # STEP-3B: NORMALIZE AND ADD TAGS TO VCF
        # ----------------------------------------------- #
        bgzip -f --threads $threads $VCFdir/${sample}_germline.vcf
        # ----------------------------------------------- #
        tabix -f --threads $threads -p vcf -0 $VCFdir/${sample}_germline.vcf.gz
         # ----------------------------------------------- #
        bcftools norm --threads $threads -c w -f ${refFasta} -m-any -o $VCFdir/${sample}_norm_germline.vcf.gz -Oz $VCFdir/${sample}_germline.vcf.gz
        # ----------------------------------------------- #
        bcftools +fill-tags --threads $threads -o $VCFdir/${sample}_norm_ft_germline.vcf.gz -Oz $VCFdir/${sample}_norm_germline.vcf.gz
        # ----------------------------------------------- #
        tabix -f --threads $threads -p vcf -0 $VCFdir/${sample}_norm_ft_germline.vcf.gz
        END_TIME "GATK_HAPLOTYPECALLER"
    else
        echo -e "${BOLD_RED}Skipping GERMLINE Variant Calling step ${NC}\n"
    fi
}   


# ----------------------------------------------- #
## STEP-4: VARIANT FILTRATION
# ----------------------------------------------- #
# HARD FILTERS FOR HAPLOTYPECALLER OUTPUT
GATK_VARIANT_FILTRATION() {
    START_TIME "GATK_VARIANT_FILTRATION"
    echo -e "${BOLD_YELLOW}>>> STEP 4 -->>> Germline Variant Filtration ${NC}\n"
    gatk VariantFiltration -R $refFasta -V $VCFdir/${sample}_norm_ft.vcf.gz -O $VCFdir/${sample}_norm_ft_filtered.vcf.gz \
        --filter-expression  'QUAL < 20.0' --filter-name "Variant_Quality_Threshold" \
        --filter-expression "MQ < 20.0" --filter-name "Mapping_Quality_Threshold" \
        --filter-expression "QD < 2.0" --filter-name "QualByDepth_Threshold" \
		--filter-expression "FS > 60.0" --filter-name "Fisher_Strand_Bias_Threshold" \
		--filter-expression "SOR > 3.0" --filter-name "Strands_Odds_Ratio_Threshold" \
		--filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum_Threshold" \
		--filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum_Threshold"
    END_TIME "GATK_VARIANT_FILTRATION"
}

# ------------------------------------------- #
## CREATE WORKFLOW AND RUN THE PIPELINE
# ------------------------------------------- #

GERMLINE() {
    NGS_ENV
    echo -e "${BOLD_BLUE}-------------------------------------------------------------------------------------------- "
    echo -e "<<<${BOLD_BLUE}      ${BOLD_PURPLE}* * * ${BOLD_YELLOW}GERMLINE VARIANT CALLING ${BOLD_PURPLE}* * *         ${BOLD_BLUE}>>>"
    echo -e "${BOLD_BLUE}-------------------------------------------------------------------------------------------- ${NC}"                                                     
    PRE_VC # Preprocessing before Variant calling
    LOGS "GATK_HAPLOTYPECALLER" # Somatic Variant calling
    LOGS "GATK_VARIANT_FILTRATION" # Somatic variant filtration
    TEMP_FILES
    # ---------------------------------- #
}

GERMLINE


# EOF
