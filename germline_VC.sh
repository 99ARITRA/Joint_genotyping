#!/usr/bin/bash
# SOF
# This pipeline is meant to do germline variant calling, filtration and annotation

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
## OUTPUR DIR PATHS
# ----------------------------------------------- #
output="/mnt/d/Bioinformatics/NGS/GENOMICS/PROJECT_OUTPUT"
bam_data=$output/BAM_DATA
prep_reports=$output/PREPROCESSING_REPORTS
vcf_data=$output/VCF_DATA
logs_dir=$output/LOG_FILES
temp=$output/TEMP # To be deleted at the end of the pipeline

# ----------------------------------------------- #
## CHECK INPUT DIR PATHS AND FILES
# ----------------------------------------------- #
inputFiles=($sample_sheet $fasta $bqsr_ref)
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
## LOG REPORTS
# ----------------------------------------------- #
LOGS() {
    processLog=$logs_dir/process.log
    $1 2>&1 | tee -a $processLog
    echo -e "\n" >> $processLog
}
# ----------------------------------------------- #   
START_TIME() {
    traceLog=$logs_dir/trace.log
    st=$(date +%H:%M:%S)
    echo "$1 START TIME: $st" >>$traceLog
}
# ----------------------------------------------- #
END_TIME() {
    et=$(date +%H:%M:%S)
    echo "$1 END TIME: $et" >>$traceLog
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

READ_ALIGNMENT() {
    echo -e "${BOLD_BLUE}>>> STEP-1 -->>> Mapping ${NC}\n"
    if [[ ! -f $bam_data/${sample}_mapped_sorted.bam || ! -f $bam_data/${sample}_bqsr.bam ]]; then
        START_TIME "ALIGNMENT"
        echo -e "${BOLD_YELLOW} Mapping to genome${NC}\n"
        bwa mem -t $threads $index $fr $rr | samtools view -1 -@ $threads -h -b -S - | \
        samtools sort -@ $threads -o $bam_data/${sample}_mapped_sorted.bam
        END_TIME "ALIGNMENT"
        # --------------------------------------------------------------- #
        echo -e "${BOLD_YELLOW} Indexing BAM${NC}\n"
        samtools index -@ $threads $bam_data/${sample}_mapped_sorted.bam
        # --------------------------------------------------------------- #
        echo -e "${BOLD_YELLOW} Generating BAM stats${NC}\n"
        samtools flagstats -@ $threads $bam_data/${sample}_mapped_sorted.bam > $prep_reports/${sample}_mapped.bam.stats
    else
        echo -e "${BOLD_RED}Skipping MAPPING step${NC}\n"
    fi
}

# ----------------------------------------------- #
## STEP-2: BAM PROCESSING
# ----------------------------------------------- #

READ_GROUPS_ADDITION() {
    echo -e "${BOLD_CYAN} Adding Read Groups${NC}\n"
    START_TIME "READ_GROUPS_ADDITION"
    gatk AddOrReplaceReadGroups -I $bam_data/${sample}_mapped_sorted.bam -O $bam_data/${sample}_rg.bam \
        --RGID rg_${sample}   --RGPL illumina  --RGSM ${sample}  --RGPU unit_${sample} --RGLB lib_${sample} -SO coordinate
    END_TIME "READ_GROUPS_ADDITION"
}
# ---------------------------------------------------------------------------------------- #
DEDUPLICATION() {
    echo -e "${BOLD_CYAN} Mark and remove Duplicates${NC}\n"
    START_TIME "DEDUPLICATION"
    gatk MarkDuplicates -I $bam_data/${sample}_rg.bam  -M $bam_data/${sample}_dup_metrics.txt  \
        -O $bam_data/${sample}_deduplicated.bam --CREATE_INDEX True
    # --------------------------------------------------------------- #
    samtools flagstats -@ $threads $bam_data/${sample}_deduplicated.bam > $prep_reports/${sample}_deduplicated.bam.stats
    END_TIME "DEDUPLICATION"
}
# ------------------------------------------------------------------------------- #
BASE_QUAL_SCORE_RECAL() {
    echo -e "${BOLD_CYAN} Base Quality Score Recalibration Step 1${NC}\n"
    START_TIME "BASE_QUAL_SCORE_RECAL"
    gatk BaseRecalibrator -I $bam_data/${sample}_deduplicated.bam \
            --known-sites $bqsr_ref \
            -O $prep_reports/${sample}_BQSR.recalibration.table \
            -R $fasta
     END_TIME "BASE_QUAL_SCORE_RECAL"
}
# ---------------------------------------------------------------------------------------- #
BQSR_APPLY() {
    echo -e "${BOLD_CYAN} Base Quality Score Recalibration Step 2${NC}\n"
    START_TIME "BQSR_APPLY"
    gatk ApplyBQSR -R $fasta -I $bam_data/${sample}_deduplicated.bam \
            --bqsr-recal-file $prep_reports/${sample}_BQSR.recalibration.table \
            -O $bam_data/${sample}_bqsr.bam -OBI
    END_TIME "BQSR_APPLY"
    # --------------------------------------------------------------- #
}


PROCESS_BAM() {
    echo -e "${BOLD_BLUE}>>> STEP-2 -->>> BAM file manipulation ${NC} \n"
    if [[ ! -f $bam_data/${sample}_deduplicated.bam  || ! -f $bam_data/${sample}_bqsr.bam ]]; then
        READ_GROUPS_ADDITION # Manipulate BAM records with sample read groups
        DEDUPLICATION # Mark and remove PCR and optical duplicates
        BASE_QUAL_SCORE_RECAL # Generate BQSR table
        BQSR_APPLY # Recalibrate Base quality scores
    else
        echo -e "${BOLD_RED}Skipping BAM PROCESSING step${NC}\n"
    fi
}


# --------------------------------------------------------------------------------------------------- #
## WORKFLOW-1: RUN ITERATION THROUGH MULTIPLE SAMPLES AND COLLECT METADATA
# -------------------------------------------------------------------------------------------------- #

NGS_PROCESSING() {
    sampleList=$(tail -n +2 $sample_sheet)

    for entry in $sampleList; do
        IFS=',' read -r SAMPLENAME R1 R2 <<< $entry
        sample=$SAMPLENAME
        fr=$R1
        rr=$R2
        
        echo -e "${BOLD_BLUE}---------------------------------------------------------------------------------------${NC}\n"
        echo -e "${BOLD_PURPLE} GERMLINE SAMPLE  :  ${BOLD_GREEN}$sample ${NC}\n"
        echo -e "${BOLD_PURPLE} REFERENCE GENOME :  ${BOLD_GREEN}$(basename $fasta .fa) ${NC}\n"
        echo -e "${BOLD_BLUE}---------------------------------------------------------------------------------------${NC}\n"

        # ---------------------------------- #
        ## PREPROCESSING PIPELINE
        # ---------------------------------- #
        BUILD_DIR
        # ----------------------------------------------- #
        LOGS "READ_ALIGNMENT"
        # ----------------------------------------------- #
        LOGS "PROCESS_BAM"
    done
}


# ----------------------------------------------- #
## STEP-3: GERMLINE VARIANT CALLING 
# ----------------------------------------------- #

GERMLINE_CALLER() {
    for bam in $(ls $bam_data/*_bqsr.bam); do
        sample=$(basename $bam _bqsr.bam)
        if [[ ! -f $vcf_data/${sample}_germline.g.vcf.gz ]]; then
            echo -e "${BOLD_YELLOW} Germline variant calling for $sample ${NC}\n"
            START_TIME "GATK_HAPLOTYPECALLER"
            gatk HaplotypeCaller -R $fasta \
                        -I $bam \
                        -O $vcf_data/${sample}_germline.g.vcf.gz \
                        -ERC GVCF
            END_TIME "GATK_HAPLOTYPECALLER"
        else
            echo -e "${BOLD_RED} Germline variant calling done for $sample ${NC}\n"
        fi
    done
}
# ----------------------------------------------- #
GVCFS_COMBINER() {
    gvcfs=( $vcf_data/*_germline.g.vcf.gz )
    # ----------------------------------------------- #
    echo -e "${BOLD_YELLOW} Combining GVCFs:\n ${BOLD_GREEN}$(ls $vcf_data/*.g.vcf.gz) ${NC} \n"
    START_TIME "COMBINE_GVCFS"
    gatk CombineGVCFs -R $fasta \
                    -V ${gvcfs[0]} \
                    -V ${gvcfs[1]} \
                    -V ${gvcfs[2]} \
                    -O $vcf_data/germline_trio.g.vcf.gz
    END_TIME "COMBINE_GVCFS"
}
# ----------------------------------------------- #
GVCF_GENOTYPER() {
    echo -e "${BOLD_YELLOW} Genotyping Trio VCF ${NC}\n"
    START_TIME "GENOTYPE_GVCFS"
    gatk GenotypeGVCFs -R $fasta \
                    -V $vcf_data/germline_trio.g.vcf.gz \
                    -O $vcf_data/joint_genotyped_trio.vcf.gz
    # ----------------------------------------------- #
    gatk IndexFeatureFile -I $vcf_data/joint_genotyped_trio.vcf.gz
    END_TIME "GENOTYPE_GVCFS"
}

# ----------------------------------------------- #
# STEP-4: NORMALIZE AND ADD TAGS TO VCF
# ----------------------------------------------- #
NORMALIZE_VCF() {
        
    if [[ ! -f $vcf_data/joint_genotyped_trio_fill-tags.vcf.gz ]]; then
        echo -e "${BOLD_YELLOW} Processing Trio VCF ${NC}\n"
        START_TIME "NORMALIZE_VCF"
        bcftools norm --threads $threads -c w -f ${fasta} -m-any -o $vcf_data/joint_genotyped_trio_normalized.vcf.gz -Oz $vcf_data/joint_genotyped_trio.vcf.gz
        # ----------------------------------------------- #
        bcftools +fill-tags --threads $threads -o $vcf_data/joint_genotyped_trio_fill-tags.vcf.gz -Oz $vcf_data/joint_genotyped_trio_normalized.vcf.gz
        # ----------------------------------------------- #
        tabix -f --threads $threads -p vcf $vcf_data/joint_genotyped_trio_fill-tags.vcf.gz
        # ----------------------------------------------- #
        echo -e "${BOLD_BLUE}---------------------------------------------------------------------------------------${NC}\n"
        count_vcf=$(bcftools view -H $vcf_data/joint_genotyped_trio_fill-tags.vcf.gz | wc -l)
        echo -e "${BOLD_GREEN}Genotyped variants: $count_vcf\n ${NC}"
        echo -e "${BOLD_BLUE}---------------------------------------------------------------------------------------${NC}\n"
        END_TIME "NORMALIZE_VCF"
    else
        echo -e "${BOLD_RED} Genotyped GVCF is normalized ${NC}\n"
    fi
}    


# --------------------------------------------------------------------------------------------------- #
## WORKFLOW-2: JOINT GENOTYOING OF TRIO-SAMPLES
# -------------------------------------------------------------------------------------------------- #

JOINT_GENOTYPING() {
    echo -e "${BOLD_BLUE}---------------------------------------------------------------------------------------${NC}\n"
    echo -e "${BOLD_PURPLE}>>> Step 3 -->>> Joint Genotyping${NC}\n"
    echo -e "${BOLD_BLUE}---------------------------------------------------------------------------------------${NC}\n"
    GERMLINE_CALLER # Call germline variants from each sample
    GVCFS_COMBINER # Combine the sample VCFs
    GVCF_GENOTYPER # Combined VCF genotyping
    NORMALIZE_VCF # VCF processing

}


# ----------------------------------------------- #
## STEP-5: VARIANT FILTRATION
# ----------------------------------------------- #
# HARD FILTERS FOR HAPLOTYPECALLER OUTPUT
HARD_FILTRATION() {
    echo -e "${BOLD_CYAN} Filtering germline variants${NC}\n"
    START_TIME "GATK_VARIANT_FILTRATION"
    gatk VariantFiltration -R $fasta -V $vcf_data/joint_genotyped_trio_fill-tags.vcf.gz -O $vcf_data/joint_genotyped_flag_filters.vcf.gz \
        --filter-expression  'FS > 60.0' --filter-name "High_Strand_Bias" \
        --filter-expression "MQ < 60.0" --filter-name "Low_Mapping_Quality" \
        --filter-expression "QD < 1.0" --filter-name "Low_Quality_by_Depth" \
        --filter-expression "MQRankSum < -10.0" --filter-name "Low_MQRankSum" \
        --filter-expression "SOR > 1.0" --filter-name "High_Odds_Ratio" \
        --filter-expression "ReadPosRankSum < -8.0" --filter-name "Low_Read_Depth" \
        --filter-expression  'QUAL < 100.0' --filter-name "Low_Quality_Variant"
    END_TIME "GATK_VARIANT_FILTRATION"
}
# ----------------------------------------------- #
SELECT_VARIANTS() {
    echo -e "${BOLD_CYAN} Eliminating variants with low genotype quality${NC}\n"
    START_TIME "SELECT_VARIANTS"
    bcftools filter -i 'FILTER="PASS"' -o $vcf_data/joint_genotyped_passed.vcf.gz -Oz $vcf_data/joint_genotyped_flag_filters.vcf.gz
    # ----------------------------------------------- #
    bcftools filter -e 'GT="./."' -o $vcf_data/joint_genotyped_filtered.vcf.gz -Oz $vcf_data/joint_genotyped_passed.vcf.gz
    # ----------------------------------------------- #
    bcftools filter -e '(GT[2]="1/1" || GT[2]="0/0")' -o $vcf_data/joint_genotyped_het.vcf.gz -Oz $vcf_data/joint_genotyped_filtered.vcf.gz
    # ----------------------------------------------- #
    bcftools filter -i 'FMT/GQ>60' -o $vcf_data/joint_genotyped_het_high_qual.vcf.gz -Oz $vcf_data/joint_genotyped_het.vcf.gz
    # ----------------------------------------------- #
    tabix -f --threads $threads -p vcf $vcf_data/joint_genotyped_het_high_qual.vcf.gz
    # ----------------------------------------------- #
    echo -e "${BOLD_BLUE}---------------------------------------------------------------------------------------${NC}\n"
    count_vcf=$(bcftools view -H $vcf_data/joint_genotyped_het_high_qual.vcf.gz | wc -l)
    echo -e "${BOLD_GREEN}Filtered variants: $count_vcf\n ${NC}"
    echo -e "${BOLD_BLUE}---------------------------------------------------------------------------------------${NC}\n"
    END_TIME "SELECT_VARIANTS"
}

# --------------------------------------------------------- #
## STEP-6: MOVE AND REMOVE THE TEMPORARY FILES
# --------------------------------------------------------- #

TEMP_FILES() {
    mv $bam_data/${sample}_mapped_sorted.bam $temp
    mv $bam_data/${sample}_mapped_sorted.bam.bai $temp
    mv $bam_data/${sample}_rg.bam $temp
    mv $bam_data/${sample}_deduplicated.bam $temp
    mv $vcf_data/germline_trio.g.vcf.gz $temp
    mv $vcf_data/joint_genotyped_trio.vcf.gz $temp
    mv $vcf_data/joint_genotyped_trio_normalized.vcf.gz $temp

    # --------------------------------------------------------- #
    rm -r $temp  
}


## --------------------------------------------------------------------------------------------------- #
## WORKFLOW-3: VARIANT FILTRATION
# -------------------------------------------------------------------------------------------------- #

GERMLINE_FILTRATION() {
    echo -e "${BOLD_BLUE}>>> STEP 4 -->>> Germline Variant Filtration ${NC}\n"
    GATK_VARIANT_FILTRATION # Applying Hard filters to germline variants
    SELECT_VARIANTS # Filtering high quality variants
    TEMP_FILES # Remove intermediate files
}


# ------------------------------------------------------------------------------ #
## MAIN WORKFLOW: EXECUTION OF WORKFLOW-1, WORKFLOW-2 AND WORKLFOW-3 FOR THE PIPELINE
# ------------------------------------------------------------------------------ #

JG_PIPELINE() {
    # NGS_ENV
    echo -e "${BOLD_BLUE}-------------------------------------------------------------------------------------------- "
    echo -e "<<<${BOLD_BLUE}      ${BOLD_PURPLE}* * * ${BOLD_YELLOW}JOINT GENOTYPING PIPELINE ${BOLD_PURPLE}* * *         ${BOLD_BLUE}>>>"
    echo -e "${BOLD_BLUE}-------------------------------------------------------------------------------------------- ${NC}"                                                     
    NGS_PROCESSING # Workflow-1: Preprocessing before Variant calling
    LOGS "JOINT_GENOTYPING" # Somatic Variant calling
    LOGS "GERMLINE_FILTRATION" # Somatic variant filtration
    echo -e "${BOLD_RED} Intermediate files have been removed ${NC} \n"
}

# ---------------------------------- #
## CLI ARGUMENTS
# --------------------------------- #

if [ $# -eq 0 ]; then
    echo -e "Please provide the following arguments
            COMMAND............
            > bash germline_VC.sh [ --samples <samplesheet.csv> ] [ --fasta <FASTA file> ] [ --index <genome index>] [ --bqsr_ref <Population VCF file>\n
                 PIPELINE PARAMETERS:
                 --samples : CSV file containing sample name,forward_fastq,reverse_fastq
                 --fasta : reference genome file in FASTA format
                 --index : genome index file created from BWA
                --bqsr_ref: Population VCF file from dbSNP"
    exit 1
fi
    
while [ $# -gt 0 ]; do
    case $1 in
    --samples)
        shift
        sample_sheet=$1 # stores sample metadata
        shift
        ;;
    --ref)
        shift
        fasta=$1 # fasta file required in analysis
        shift
        ;;
    --idx)
        shift
        index=$1 # genome index file requried in mapping step
        shift
        ;;
    --bqsr_ref)
        shift
        bqsr_ref=$1 # Population VCF file requried in BQSR step
        shift
        ;;
    --cpus)
        shift
        threads=$1
        shift
        ;;
    * | -h)
        echo -e "Wrong argument entered........\n
                COMMAND............
                > bash germline_VC.sh [ --samples <samplesheet.csv> ] [ --fasta <FASTA file> ] [ --index <genome index>] [ --bqsr_ref <Population VCF file>\n
                    PIPELINE PARAMETERS:
                    --samples : CSV file containing sample name,forward_fastq,reverse_fastq
                    --ref : reference genome file in FASTA format
                    --idx : genome index file created from BWA
                    --bqsr_ref : Population VCF file from dbSNP
                    --threads : No. of CPUs to provide in process"
        exit 1
    esac
done

# ---------------------------------- #
## INITIATE PIPELINE RUN
# --------------------------------- #
JG_PIPELINE

# EOF