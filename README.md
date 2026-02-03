# Joint genotyping
This is a project on joint genotyping using NGS tools on a father, mother and proband sample for identification of variants related to Osteopertrosis. The analysis method was adapted from a Galaxy tutorial project on the disease.
## Tools used
- BWA-MEM
- SAMTOOLS
- GATK BAM PROCESSING TOOLS
  - ADDORREPLACEREADGROUPS
  - MARKDUPLICATES
  - BASEQUALITYSCORERECALIBRATION
  - APPLYBQSR
- VARIANT CALLING
  - HAPLOTYPECALLER (GVCF MODE)
  - COMBINEGVCFS
  - GENOTYPEGVCFS
- VCF PROCESSING
  - BCFTOOLS NORM
  - BCFTOOLS FILLTAGS
- VARIANT FILTRATION
  - GATK VARIANTFILTRATION
 
## Installation instructions
    1. Download the conda executable < wget https://github.com/conda-forge/miniforge/releases/download/25.11.0-1/Miniforge3-Linux-x86_64.sh >
    2. Change permissions < chmod 777 Miniforge3-Linux-x86_64.sh >
    3. Run < ./Miniforge3-Linux-x86_64.sh >
    4. Set up the conda environment.
    5. Create a new conda environment < conda create -n < envname > >
    6. Install the tools < mamba install bwa-mem samtools gatk4 bcftools >
    7. Download the repo < git clone <repo_name> >
    8. Navigate to the directory in terminal and run the script as shown below.
    
## Pipeline execution
 To execute the bash script type the following commands in a linux terminal or WSL terminal:
 `bash germline_VC.sh`
