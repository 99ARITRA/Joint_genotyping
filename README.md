# Joint genotyping

This NGS data analysis pipeline is based on joint genotyping.  The analysis method was adapted from the  Galaxy tutorial project <https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/exome-seq/tutorial.html>. The project involves identification of variants in a proband sample and its parents and comparing their genotypes to predict the probable causative variants for the genetic disorder in the proband. The pipeline is designed to be run on a Linux terminal or WSL terminal. The pipeline takes as input a samplesheet file in CSV format. The pipeline outputs BAM, VCF and a table of variants with their annotations.

## Tools used

- BWA-MEM 2 (MAPPING)
- SAMBAMBA (BAM FILE PROCESSING)
- FREEBAYES ( JOINT VARIANT CALLING)
- BCFTOOLS (VCF FILE PROCESSING)
- SNPEFF (VARIANT ANNOTATION)
- VCFANNO (VARIANT ANNOTATION)
- TABIX

## Installation instructions

  1. Download Pixi: ` curl -fsSL https://pixi.sh/install.sh | sh `
  2. Add the following to ~/.bashrc: `eval "$(pixi completion --shell bash)"` to enable autocompletion. Restart the shell.
  3. Install the tools using Pixi: `pixi global install [Tool name]`. Refer the section [Tools used](#tools-used)
  4. Download tabix >> ` sudo apt install tabix `
  5. Install pandas library >> ` sudo apt install python3-pandas OR pip install pandas `.

  6. Download the repo: ` git clone https://github.com/99ARITRA/Joint_genotyping.git `
  7. Navigate to the directory "Joint_genotyping" in terminal and run the script as shown below.

## Pipeline execution

 To execute the bash script type the following commands in a linux terminal or WSL terminal:
 `bash joint_genotyping.sh --help`. It will display the list of arguments to be supplied to the command and how to run the command in a linux terminal.

## Instructions

  The samplesheet layout of <samplesheet.csv> for the pipeline is as follows:

  | COHORT_NAME | SAMPLE_NAME | R1 | R2 |
  | ------------- | ------------- | ---- | ---- |
  | ABC | proband | file_R1.fastq.gz | file_R2.fastq.gz |
  | ABC | parent1 | file_R1.fastq.gz | file_R2.fastq.gz |
  | ABC | parent2 | file_R1.fastq.gz | file_R2.fastq.gz |

  For joint genotyping, three samples need to used for the analysis: one for the proband sample, and the other two for each of the parent.
