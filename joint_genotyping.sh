#!/usr/bin/bash
# SOF
# This pipeline is meant to do germline variant calling, filtration and annotation

set -euo pipefail # Pipeline falure command

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
## CREATE THE OUTPUT DIRECTORIES
# ----------------------------------------------- #

BUILD_DIR() {
    output=/output/"JG_OUTPUT"/$cohort
    bam_data=$output/BAM_DATA
    prep_reports=$output/PREPROCESSING_REPORTS
    vcf_data=$output/VCF_DATA
    logs_dir=$output/LOG_FILES
    temp=$output/TEMP # To be deleted at the end of the pipeline
    mkdir -p $output # output reports and files directory
    mkdir -p $bam_data # Storing BAM files
    mkdir -p $prep_reports # Storing preprocessing reports
    mkdir -p $vcf_data # Storing unfiltered and filtered VCF files
    mkdir -p $logs_dir # Log report
    mkdir -p $temp # Directory to store intermediate files in pipeline. Later to be deleted during pipeline run
    # ----------------------------------------------------- #
    processLog=$logs_dir/process.log # Log report file
    echo -e "\n============================= LOG REPORT =============================\n" > "$processLog"
}


# ----------------------------------------------- #
## STEP-1: MAPPING
# ----------------------------------------------- #

READ_ALIGNMENT() {
    echo -e "${BOLD_BLUE}>>> STEP-1 -->>> Mapping | SAMPLE: ${BOLD_PURPLE}$sample ${NC}\n"
    bamFile=$bam_data/${sample}_mapped_sorted.bam # mapped bam file
    if [[ ! -f $bamFile ]]; then
        echo -e "${BOLD_YELLOW} Mapping to : ${BOLD_GREEN}$(basename $fasta .fa) ${NC}\n"
        bwa-mem2 mem -v 1 -R $read_group -t $threads $index $R1 $R2 | sambamba view -S -f bam -l 0 -t $threads -p /dev/stdin | \
        sambamba sort -m 450MB -t $threads -l 0 -p -o $bamFile /dev/stdin
        # --------------------------------------------------------------- #
        echo -e "${BOLD_YELLOW} Indexing BAM${NC}\n"
        sambamba index -t $threads -p $bamFile
        # --------------------------------------------------------------- #
        echo -e "${BOLD_YELLOW} Generating BAM stats${NC}\n"
        sambamba flagstat -t $threads -p $bamFile > $prep_reports/${sample}_mapped.bam.stats
    else
        echo -e "${BOLD_RED}Skipping MAPPING step${NC}\n"
    fi
}

# ----------------------------------------------- #
## STEP-2: BAM PROCESSING
# ----------------------------------------------- #

DEDUPLICATION() {
    echo -e "${BOLD_BLUE}>>> STEP-2 -->>> Deduplication | SAMPLE: ${BOLD_PURPLE}$sample ${NC}\n"
    dedupBam=$bam_data/${sample}_deduplicated_sorted.bam # deduplicated bam file
    if [[ ! -f $dedupBam ]]; then
        echo -e "${BOLD_YELLOW} Marking and removing duplicates from BAM${NC}\n"
        sambamba markdup -r -t $threads -p  $bamFile /dev/stdout | \
        sambamba sort -m 450MB -t $threads -l 9 -p -o $dedupBam /dev/stdin
        # --------------------------------------------------------------- #
        echo -e "${BOLD_YELLOW} Indexing BAM${NC}\n"
        sambamba index -t $threads -p $dedupBam
        # --------------------------------------------------------------- #
        echo -e "${BOLD_YELLOW} Generating BAM stats${NC}\n"
        sambamba flagstat -t $threads -p $dedupBam > $prep_reports/${sample}_deduplicated.bam.stats
    else
        echo -e "${BOLD_RED} Skipping DEDUPLICATION step${NC}\n"
    fi
}

# ----------------------------------------------- 3#
## STEP-3: TEXT FILE WITH PATH TO TRIO BAM FILES
# ----------------------------------------------- #

TRIPLE_BAM() {
    find $bam_data -type f -name *_deduplicated_sorted.bam -exec realpath {} \; > $bam_data/trio_bams.txt
}


# ------------------------------------------------------------------------------------------------------------- #
## WORKFLOW-1: RUN ITERATION THROUGH MULTIPLE SAMPLES AND COLLECT METADATA (ALIGNMENT + BAM PROCESSING)
# ------------------------------------------------------------------------------------------------------------ #

NGS_PROCESSING() {       
    echo -e "${BOLD_CYAN}===============================================================================================================================================================================================${NC}"
    echo -e "${BOLD_PURPLE}SAMPLE METADATA:"
    echo -e "${BOLD_GREEN}$(column -t -s "," $sample_sheet)"
    echo -e "${BOLD_CYAN}===============================================================================================================================================================================================${NC}\n"    
    
    sampleList=$(tail -n +2 $sample_sheet)

    for entry in $sampleList; do
        IFS=',' read -r COHORT_NAME SAMPLE_NAME R1 R2 <<< $entry
        cohort=$COHORT_NAME
        sample=$SAMPLE_NAME
        fr=$R1
        rr=$R2
        read_group="@RG\tID:$sample\tSM:$sample\tPL:ILLUMINA"
        # ---------------------------------- #
        ## PREPROCESSING PIPELINE
        # ---------------------------------- #
        BUILD_DIR # Create directories for storing results         
        READ_ALIGNMENT 2>> $processLog # Mapping of reads to reference genome
        DEDUPLICATION 2>> $processLog # Marking and removing duplicates from BAM
    done
    # ----------------------------------------------- #
    TRIPLE_BAM # Create a text file with path to trio bam files
}


# ----------------------------------------------- #
## STEP-4: GERMLINE VARIANT CALLING 
# ----------------------------------------------- #

GERMLINE_CALLER() {
    echo -e "${BOLD_YELLOW} Germline variant calling | BAM files:\n${BOLD_GREEN}$(cat $bam_data/trio_bams.txt) ${NC}"
    echo -e "${BOLD_CYAN}===============================================================================================================================================================================================${NC}\n"
    freebayes -f $fasta -0 --genotype-qualities -L $bam_data/trio_bams.txt > $trioVcf
    # ----------------------------------------------- #
    bgzip -f -@ $threads $trioVcf
    # ----------------------------------------------- #
    tabix -f -p vcf $trioVcf.gz
}

# ----------------------------------------------- #
# STEP-5: NORMALIZE AND ADD TAGS TO VCF
# ----------------------------------------------- #

NORMALIZE_VCF() {
    echo -e "${BOLD_YELLOW} Processing Trio VCF ${NC}\n"
    bcftools norm --threads $threads -c w -f ${fasta} -m-any -o $normVcf -Oz $trioVcf.gz
    # ----------------------------------------------- #
    bcftools +fill-tags --threads $threads -o $ftVcf -Oz $normVcf
    # ----------------------------------------------- #
    tabix -f -p vcf $ftVcf
    # ----------------------------------------------- #
    echo -e "${BOLD_CYAN}===============================================================================================================================================================================================${NC}\n"
    count_vcf=$(bcftools view -H $ftVcf | wc -l)
    echo -e "${BOLD_GREEN} Variants called: $count_vcf\n ${NC}"
    echo -e "${BOLD_CYAN}===============================================================================================================================================================================================${NC}\n"
}    

# ----------------------------------------------- #
## STEP-6: VARIANT FILTRATION
# ----------------------------------------------- #

FILTER_VCF() {
    echo -e "${BOLD_YELLOW} Filtering variants based on ${BOLD_PURPLE}Hard Filters${NC}\n"
    bcftools filter --threads $threads -i '(QUAL>30) && INFO/DP>30' -o $filtVcf1 \
                        -Oz $ftVcf
    # ----------------------------------------------- #
    echo -e "${BOLD_YELLOW} Filtering variants based on ${BOLD_PURPLE}Genotype${NC}\n"
    bcftools filter --threads $threads -e 'FORMAT/GT[2]="mis" || FORMAT/GT[0]="mis" || FORMAT/GT[1]="mis"' -o $filtVcf2 \
                        -Oz $filtVcf1
    # ----------------------------------------------- #
    tabix -f -p vcf $filtVcf2
    # ----------------------------------------------- #
    echo -e "${BOLD_YELLOW} Generating VCF stats ${NC}\n"
    bcftools stats -s - --threads $threads $filtVcf2 > $prep_reports/$(basename $filtVcf2 .gz).stats
}


# --------------------------------------------------------------------------------------------------- #
## WORKFLOW-2: JOINT GENOTYOING OF TRIO-SAMPLES
# -------------------------------------------------------------------------------------------------- #

JOINT_GENOTYPING() {
    trioVcf=$vcf_data/trio_jg.vcf # joint genotyped trio vcf file
    normVcf=$vcf_data/trio_jg_normalized.vcf.gz # normalized vcf file
    ftVcf=$vcf_data/trio_jg_fill-tags.vcf.gz # fill-tags vcf file
    filtVcf1=$vcf_data/jg_passed_v1.vcf.gz # filtered vcf file 1
    filtVcf2=$vcf_data/jg_passed_v2.vcf.gz # filtered vcf file 2
    echo -e "${BOLD_BLUE}>>> Step 3 -->>> Joint Genotyping${NC}\n"
    if [[ ! -f $ftVcf || ! -f $filtVcf2 ]]; then
        GERMLINE_CALLER 2>> $processLog # Call germline variants from each sample
        NORMALIZE_VCF 2>> $processLog # VCF processing
        FILTER_VCF 2>> $processLog # Filtering high quality variants
    else
        echo -e "${BOLD_RED} Genotyped GVCF is normalized ${NC}\n"
    fi
}


# ----------------------------------------------- #
## STEP-7: VARIANT ANNOTATION
# ----------------------------------------------- #

SNPEFF() {
    annotVcf1=$vcf_data/jg_annotated_v1.vcf # annotated vcf file 1
    annotVcf2=$vcf_data/jg_annotated_v2.vcf # annotated vcf file 2
    annotVcf3=$vcf_data/jg_annotated_v3.vcf # annotated vcf file 3
    annotVcf4=$vcf_data/jg_annotated_v4.vcf # annotated vcf file 4
    echo -e "${BOLD_YELLOW} Annotating VCF | ${BOLD_PURPLE}SnpEff ${NC}\n"
    snpEff ann -Xmx${memo}g -i vcf -o vcf  -htmlStats $vcf_data/annotation.html \
                        -hgvs -canon -no-downstream -no-intergenic -no-intron -no-upstream -no-utr \
                        $ann_genome $filtVcf2 > $annotVcf1
    # --------------------------------------------------------------------------------- #
    SnpSift Annotate -id $dbsnp $annotVcf1 > $annotVcf2
    # --------------------------------------------------------------------------------- #
    SnpSift Filter "(ANN[*].IMPACT = 'HIGH') | (ANN[*].IMPACT = 'MODERATE')" $annotVcf2 > $annotVcf3
}
# ---------------------------------------------------- #
CLINVAR() {
    clinvarLua="clinvar.lua"
    clinvarToml="clinvar.toml"
    echo -e "${BOLD_YELLOW} Annotating VCF |${BOLD_PURPLE} ClinVar${NC}\n"
    vcfanno -p $threads -lua $clinvarLua $clinvarToml \
                                $annotVcf3 > $annotVcf4
}
# ---------------------------------------------------- #
TABULATION() {
    tab1=$output/High-confidence_variants.tsv # table 1
    tab2=$output/genotyped_variants.tsv # table 2
    tab3=$output/Final_variants_table.tsv # table 3
    echo -e "${BOLD_YELLOW} Collecting SNPEFF fields ${NC}\n"
    SnpSift extractFields $annotVcf4 \
                        CHROM POS ID REF ALT \
                        ANN[*].EFFECT ANN[*].IMPACT ANN[*].GENE ANN[*].BIOTYPE ANN[*].FEATUREID \
                        ANN[*].HGVS_C ANN[*].HGVS_P \
                        AF_EXAC CLNSIG CLNHGVS > $tab1
    # ---------------------------------------------------- #
    echo -e "${BOLD_YELLOW} Collecting Genotype Fields ${NC}\n"
    gatk --java-options "-Xmx${memo}g" VariantsToTable -V $annotVcf4 \
                                                    -O $tab2 \
                                                    -F CHROM -F POS -F ID -F REF -F ALT -GF GT
    # ---------------------------------------------------- #
    echo -e "${BOLD_YELLOW} Merging tables${NC}\n"
    python merge_tables.py $tab1 $tab2 $tab3
}


# --------------------------------------------------------- #
## REMOVE THE TEMPORARY FILES
# --------------------------------------------------------- #

TEMP_FILES() {
    mv $normVcf $temp
    mv $filtVcf1 $temp
    mv $annotVcf1 $temp
    mv $annotVcf2 $temp
    mv $annotVcf3 $temp
    mv $tab1 $temp
    mv $tab2 $temp
    # --------------------------------------------------------- #
    rm -r $temp
    echo -e "${BOLD_RED} Intermediate files have been removed ${NC} \n"
}


## --------------------------------------------------------------------------------------------------- #
## WORKFLOW-3: VARIANT ANNOTATION
# -------------------------------------------------------------------------------------------------- #

VARIANT_ANNOTATION() {
    echo -e "${BOLD_BLUE}>>> STEP 4 -->>> Variant Annotation ${NC}\n"
    SNPEFF 2>> $processLog # Annotate VCF with SnpEff
    CLINVAR 2>> $processLog # Annotate VCF with ClinVar database information
    TABULATION 2>> $processLog # Create table from VCF file
    TEMP_FILES 2>> $processLog # Remove intermediate files
}


# ------------------------------------------------------------------------------ #
## MAIN WORKFLOW: EXECUTION OF WORKFLOW-1, WORKFLOW-2 AND WORKLFOW-3 FOR THE PIPELINE
# ------------------------------------------------------------------------------ #

JG_PIPELINE() {
    echo -e "${BOLD_BLUE} =========================================================================================================================== "
    echo -e "   <<<${BOLD_BLUE}      ${BOLD_PURPLE}* * * ${BOLD_YELLOW}JOINT GENOTYPING PIPELINE ${BOLD_PURPLE}* * *         ${BOLD_BLUE}>>>"
    echo -e "${BOLD_BLUE} =========================================================================================================================== ${NC}"
    # ----------------------------------------------- #
    # MAIN WORKFLOW
    # ----------------------------------------------- #                                              
    NGS_PROCESSING # Workflow-1: Preprocessing of NGS data
    JOINT_GENOTYPING # Workflow-2: Joint Genotyping of three samples
    VARIANT_ANNOTATION # Workflow-3: Germline Variant Annotation
}

# ---------------------------------- #
## CLI ARGUMENTS
# --------------------------------- #

if [ $# -eq 0 ]; then
    echo -e "Please provide the following arguments
            COMMAND............
            > bash joint_genotyping.sh [ --samples <samplesheet.csv> ] [ --ref <FASTA file> ] [ --idx <genome index> ] [ --annGenome <genome version for variant annotation> ]  [ --dbsnp <Population VCF file> ] [ --clinvar <ClinVar VCF file> ] [ --cpus <cpus> ] [ --mem <memory in GB> ]\n
                 PIPELINE PARAMETERS:
                    --samples : CSV file containing sample name,forward_fastq,reverse_fastq
                    --ref : reference genome file in FASTA format
                    --idx : genome index file created from BWA-MEM2
                    --annGenome : Specify the genome version for SnpEff annotation
                    --dbsnp : Population VCF file from dbSNP
                    --clinvar : ClinVar VCF file required in annotation step
                    --cpus : No. of CPUs to provide in process
                    --mem : Memory allocation for GATK tools"
    exit 1
fi
    
while [ $# -gt 0 ]; do
    case $1 in
    --samples)
        shift
        sample_sheet=$1 # stores sample metadata [CSV format]
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
    --annGenome)
        shift
        ann_genome=$1 # Enable BQSR step in pipeline [value: yes/ no]
        shift
        ;;
    --dbsnp)
        shift
        dbsnp=$1 # Population VCF file requried in BQSR step
        shift
        ;;
    --clinvar)
        shift
        clinvar=$1 # ClinVar VCF file required in annotation step
        shift
        ;;
    --cpus)
        shift
        threads=$1 # No. of CPUs required in pipeline
        shift
        ;;
    --mem)
        shift
        memo=$1 # Memory allocated for GATK Tools
        shift
        ;;
    --help | -h)
        echo -e "Wrong argument entered........\n
                COMMAND............
                > bash joint_genotyping.sh [ --samples <samplesheet.csv> ] [ --ref <FASTA file> ] [ --idx <genome index> ] [--annGenome <genome version for variant annotation>] [ --dbsnp <Population VCF file> ] [--clinvar <ClinVar VCF file>] [ --cpus <cpus> ] [ --mem <memory in GB> ]\n
                    PIPELINE PARAMETERS:
                    --samples : CSV file containing sample name,forward_fastq,reverse_fastq
                    --ref : reference genome file in FASTA format
                    --idx : genome index file created from BWA
                    --annGenome : Specify the genome version for SnpEff annotation
                    --dbsnp : Population VCF file from dbSNP
                    --clinvar : ClinVar VCF file required in annotation step
                    --cpus : No. of CPUs to provide in process
                    --mem : Memory allocation for GATK tools (in GB)"
        exit 1
    esac
done

# ---------------------------------- #
## INITIATE PIPELINE RUN
# --------------------------------- #
JG_PIPELINE

# EOF