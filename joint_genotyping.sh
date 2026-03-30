#!/usr/bin/bash
# SOF
# This pipeline is meant to do germline variant calling, filtration and annotation

# ----------------------------------------------- #
## CLI COLOURS 
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
bam_data=$output/BAM_DATA
prep_reports=$output/PREPROCESSING_REPORTS
vcf_data=$output/VCF_DATA
logs_dir=$output/LOG_FILES
temp=$output/TEMP # To be deleted at the end of the pipeline

# ----------------------------------------------- #
## PIPELINE ENVIRONNMENT ACTIVATION
# ----------------------------------------------- #
CONDA_ACTIVATION() {
    source $HOME/miniforge3/etc/profile.d/conda.sh
    conda activate $conda_env
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
    echo -e "${BOLD_BLUE}>>> STEP-1 -->>> Mapping | SAMPLE: ${BOLD_PURPLE}$sample ${NC}\n"
    if [[ ! -f $bam_data/${sample}_mapped_sorted.bam || ! -f $bam_data/${sample}_bqsr.bam ]]; then
        echo -e "${BOLD_YELLOW} Mapping to genome ${BOLD_GREEN}$(basename $fasta .fa) ${NC}\n"
        bwa-mem2 mem -v 1 -t $threads $index $fr $rr | samtools view -1 -@ $threads -h -b -S - | \
        samtools sort -@ $threads -o $bam_data/${sample}_mapped_sorted.bam
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
    gatk --java-options "-Xmx${memo}g" AddOrReplaceReadGroups -I $bam_data/${sample}_mapped_sorted.bam -O $bam_data/${sample}_rg.bam \
        --RGID rg_${sample}   --RGPL illumina  --RGSM ${sample}  --RGPU unit_${sample} --RGLB lib_${sample} -SO coordinate --VERBOSITY ERROR
}
# ---------------------------------------------------------------------------------------- #
DEDUPLICATION() {
    echo -e "${BOLD_CYAN} Mark and remove Duplicates${NC}\n"
    sambamba-markdup -r -t $threads -p $bam_data/${sample}_rg.bam $bam_data/${sample}_deduplicated.bam
    # --------------------------------------------------------------- #
    sambamba-sort -m 450MB -t $threads -p -o $bam_data/${sample}_deduplicated_sorted.bam $bam_data/${sample}_deduplicated.bam
    # --------------------------------------------------------------- #    
    sambamba-index -t $threads -p $bam_data/${sample}_deduplicated_sorted.bam  $bam_data/${sample}_deduplicated_sorted.bam.bai
    # --------------------------------------------------------------- #
    sambamba-flagstat -t $threads -p $bam_data/${sample}_deduplicated_sorted.bam > $prep_reports/${sample}_deduplicated.bam.stats
}
# ------------------------------------------------------------------------------- #
BASE_QUAL_SCORE_RECAL_1() {
    echo -e "${BOLD_CYAN} Base Quality Score Recalibration Step 1${NC}\n"
    gatk --java-options "-Xmx${memo}g" BaseRecalibrator -I $bam_data/${sample}_deduplicated_sorted.bam \
            --known-sites $bqsr_ref \
            -O $prep_reports/${sample}_BQSR.recalibration.table \
            -R $fasta --verbosity ERROR
}
# ---------------------------------------------------------------------------------------- #
BASE_QUAL_SCORE_RECAL_2() {
    echo -e "${BOLD_CYAN} Base Quality Score Recalibration Step 2${NC}\n"
    gatk --java-options "-Xmx${memo}g" ApplyBQSR -R $fasta -I $bam_data/${sample}_deduplicated_sorted.bam \
            --bqsr-recal-file $prep_reports/${sample}_BQSR.recalibration.table \
            -O $bam_data/${sample}_bqsr.bam -OBI --verbosity ERROR
    # --------------------------------------------------------------- #
    samtools flagstats -@ $bam_data/${sample}_bqsr.bam > $prep_reports/${sample}_bqsr.bam.stats
}


PROCESS_BAM() {
    echo -e "${BOLD_BLUE}>>> STEP-2 -->>> BAM file manipulation  | SAMPLE: ${BOLD_PURPLE}$sample ${NC} \n"
    if [[ ! -f $bam_data/${sample}_deduplicated.bam  || ! -f $bam_data/${sample}_bqsr.bam ]]; then
        READ_GROUPS_ADDITION # Manipulate BAM records with sample read groups
        DEDUPLICATION # Mark and remove PCR and optical duplicates
        BASE_QUAL_SCORE_RECAL_1 # Generate BQSR table
        BASE_QUAL_SCORE_RECAL_2 # Recalibrate Base quality scores
    else
        echo -e "${BOLD_RED}Skipping BAM PROCESSING step${NC}\n"
    fi
}


# ------------------------------------------------------------------------------------------------------------- #
## WORKFLOW-1: RUN ITERATION THROUGH MULTIPLE SAMPLES AND COLLECT METADATA (ALIGNMENT + BAM PROCESSING)
# ------------------------------------------------------------------------------------------------------------ #

NGS_PROCESSING() {       
    echo -e "${BOLD_BLUE}---------------------------------------------------------------------------------------${NC}"
    echo -e "${BOLD_PURPLE}SAMPLE METADATA:"
    echo -e "${BOLD_GREEN}$(column -t -s "," $sample_sheet)"
    echo -e "${BOLD_BLUE}---------------------------------------------------------------------------------------${NC}\n"    
    
    sampleList=$(tail -n +2 $sample_sheet)

    for entry in $sampleList; do
        IFS=',' read -r SAMPLENAME R1 R2 <<< $entry
        sample=$SAMPLENAME
        fr=$R1
        rr=$R2
        # ---------------------------------- #
        ## PREPROCESSING PIPELINE
        # ---------------------------------- #
        READ_ALIGNMENT
        # ----------------------------------------------- #
        PROCESS_BAM
    done
}


# ----------------------------------------------- #
## STEP-3: GERMLINE VARIANT CALLING 
# ----------------------------------------------- #

GERMLINE_CALLER() {
    for bam in $(ls $bam_data/*_bqsr.bam); do
        sample=$(basename $bam _bqsr.bam)
        if [[ ! -f $vcf_data/${sample}_germline.g.vcf.gz ]]; then
            echo -e "${BOLD_YELLOW} Germline variant calling | SAMPLE: ${BOLD_PURPLE}$sample ${NC}\n"
            gatk --java-options "-Xmx${memo}g" HaplotypeCaller -R $fasta \
                        -I $bam \
                        -O $vcf_data/${sample}_germline.g.vcf.gz \
                        -ERC GVCF --verbosity ERROR
        else
            echo -e "${BOLD_RED} Germline variant calling done for $sample ${NC}\n"
        fi
    done
}
# ----------------------------------------------- #
GVCFS_COMBINER() {
    gvcfs=( $vcf_data/*_germline.g.vcf.gz )
    # ----------------------------------------------- #
    echo -e "${BOLD_YELLOW} Combining GVCFs:\n ${BOLD_GREEN}$(ls $vcf_data/*_germline.g.vcf.gz) ${NC} \n"
    gatk --java-options "-Xmx${memo}g" CombineGVCFs -R $fasta \
                    -V ${gvcfs[0]} \
                    -V ${gvcfs[1]} \
                    -V ${gvcfs[2]} \
                    -O $vcf_data/germline_trio.g.vcf.gz --verbosity ERROR
}
# ----------------------------------------------- #
GVCF_GENOTYPER() {
    echo -e "${BOLD_YELLOW} Genotyping Trio VCF ${NC}\n"
    gatk --java-options "-Xmx${memo}g" GenotypeGVCFs -R $fasta \
                    -V $vcf_data/germline_trio.g.vcf.gz \
                    -O $vcf_data/joint_genotyped_trio.vcf.gz --verbosity ERROR
    # ----------------------------------------------- #
    gatk --java-options "-Xmx${memo}g" IndexFeatureFile -I $vcf_data/joint_genotyped_trio.vcf.gz
}

# ----------------------------------------------- #
# STEP-4: NORMALIZE AND ADD TAGS TO VCF
# ----------------------------------------------- #

NORMALIZE_VCF() {
    if [[ ! -f $vcf_data/joint_genotyped_trio_fill-tags.vcf.gz ]]; then
        echo -e "${BOLD_YELLOW} Processing Trio VCF ${NC}\n"
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
FILTER_VCF() {
    echo -e "${BOLD_CYAN} Filtering out germline variants with low depth, low genotype quality and low quality-by-depth${NC}\n"
    gatk --java-options "-Xmx${memo}g" FilterVcf -R $fasta -I $vcf_data/joint_genotyped_trio_fill-tags.vcf.gz \
                                            -O $vcf_data/joint_genotyped_flagged.vcf.gz \
                                            --MIN_DP 30.0 \
                                            --MIN_QD 1.0 \
                                            --MIN_GQ 60.0 \
                                            --VERBOSITY ERROR
}
# ----------------------------------------------- #
SELECT_VARIANTS() {
    echo -e "${BOLD_CYAN} Eliminating variants with low variant quality, mapping quality${NC}\n"
    bcftools filter -i 'FILTER="PASS" && QUAL>30 && INFO/MQ>60' -o $vcf_data/joint_genotyped_first_pass.vcf.gz -Oz $vcf_data/joint_genotyped_flagged.vcf.gz
    # ----------------------------------------------- #
    bcftools filter -e 'FMT/GT[0]=="mis" && FMT/GT[1]=="mis" && FMT/GT[2]=="mis"' -o $vcf_data/joint_genotyped_second_pass.vcf.gz -Oz $vcf_data/joint_genotyped_first_pass.vcf.gz
    # ----------------------------------------------- #
    tabix -f --threads $threads -p vcf $vcf_data/joint_genotyped_second_pass.vcf.gz
    # ----------------------------------------------- #
    echo -e "${BOLD_BLUE}---------------------------------------------------------------------------------------${NC}\n"
    count_vcf=$(bcftools view -H $vcf_data/joint_genotyped_second_pass.vcf.gz | wc -l)
    echo -e "${BOLD_GREEN}Filtered variants: $count_vcf\n ${NC}"
    echo -e "${BOLD_BLUE}---------------------------------------------------------------------------------------${NC}\n"
}

# ----------------------------------------------- #
## STEP-6: VARIANT FILTRATION
# ----------------------------------------------- #
# CONVERTING VCF TO TSV FORMAT
TABULATION() {
    echo -e "${BOLD_CYAN} Tabulating VCF${NC}\n"
    gatk --java-options "-Xmx${memo}g" VariantsToTable -V $vcf_data/joint_genotyped_second_pass.vcf.gz \
                                                    -O $vcf_data/High-confidence_variants.tsv \
                                                    -F CHROM -F POS -F REF -F ALT -GF GT
}

# --------------------------------------------------------- #
## STEP-6: MOVE AND REMOVE THE TEMPORARY FILES
# --------------------------------------------------------- #

TEMP_FILES() {
    mv $bam_data/*_mapped_sorted.bam* $temp
    mv $bam_data/*_rg.bam $temp
    mv $bam_data/*_deduplicated.bam $temp
    mv $bam_data/*_deduplicated_sorted.bam* $temp
    mv $vcf_data/*_germline.g.vcf.gz* $temp
    mv $vcf_data/germline_trio.g.vcf.gz* $temp
    mv $vcf_data/joint_genotyped_trio.vcf.gz* $temp
    mv $vcf_data/joint_genotyped_trio_normalized.vcf.gz $temp
    mv $vcf_data/joint_genotyped_trio_fill-tags.vcf.gz $temp
    mv $vcf_data/joint_genotyped_first_pass.vcf.gz $temp      
    # --------------------------------------------------------- #
    rm -r $temp  
}


## --------------------------------------------------------------------------------------------------- #
## WORKFLOW-3: VARIANT FILTRATION
# -------------------------------------------------------------------------------------------------- #

GERMLINE_FILTRATION() {
    echo -e "${BOLD_BLUE}>>> STEP 4 -->>> Germline Variant Filtration ${NC}\n"
    FILTER_VCF # Applying Hard filters to germline variants
    SELECT_VARIANTS # Filtering high quality variants
    TEMP_FILES # Remove intermediate files
}


# ------------------------------------------------------------------------------ #
## MAIN WORKFLOW: EXECUTION OF WORKFLOW-1, WORKFLOW-2 AND WORKLFOW-3 FOR THE PIPELINE
# ------------------------------------------------------------------------------ #

JG_PIPELINE() {
    echo -e "${BOLD_BLUE}-------------------------------------------------------------------------------------------- "
    echo -e "<<<${BOLD_BLUE}      ${BOLD_PURPLE}* * * ${BOLD_YELLOW}JOINT GENOTYPING PIPELINE ${BOLD_PURPLE}* * *         ${BOLD_BLUE}>>>"
    echo -e "${BOLD_BLUE}-------------------------------------------------------------------------------------------- ${NC}"                                                     
    LOGS "NGS_PROCESSING" # Workflow-1: Preprocessing of NGS data
    LOGS "JOINT_GENOTYPING" # Workflow-2: Joint Genotyping of three samples
    LOGS "GERMLINE_FILTRATION" # Workflow-3: Germline Variant Filtration
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
                    --ref : reference genome file in FASTA format
                    --idx : genome index file created from BWA-MEM2
                    --bqsr_ref : Population VCF file from dbSNP
                    --cpus : No. of CPUs to provide in process
                    --gatk_mem : Memory allocation for GATK tools
                    --conda_env : Initiate the conda environment containing the tools
                    --output : Output directory to store results"
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
        threads=$1 # No. of CPUs required in pipeline
        shift
        ;;
    --gatk_mem)
        shift
        memo=$1 # Memory allocated for GATK Tools
        shift
        ;;
    --conda_env)
        shift
        conda_env=$1
        CONDA_ACTIVATION
        shift
        ;;
    --output)
        shift
        output=$1
        BUILD_DIR
        shift
        ;;
    * | -h)
        echo -e "Wrong argument entered........\n
                COMMAND............
                > bash germline_VC.sh [ --samples <samplesheet.csv> ] [ --ref <FASTA file> ] [ --idx <genome index>] [ --bqsr_ref <Population VCF file> ]  [ --cpus <cpus> ] [ --gatk_mem <memory in GB> ] [ --conda_env <env_name> ] [ --output <output directory>\n
                    PIPELINE PARAMETERS:
                    --samples : CSV file containing sample name,forward_fastq,reverse_fastq
                    --ref : reference genome file in FASTA format
                    --idx : genome index file created from BWA
                    --bqsr_ref : Population VCF file from dbSNP
                    --cpus : No. of CPUs to provide in process
                    --gatk_mem : Memory allocation for GATK tools
                    --conda_env : Activate the conda environment containing the tools
                    --output : Output directory to store results"

        exit 1
    esac
done

# ---------------------------------- #
## INITIATE PIPELINE RUN
# --------------------------------- #
JG_PIPELINE

# EOF