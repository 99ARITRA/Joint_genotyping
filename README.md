# Joint genotyping

This is a project on joint genotyping of a father, mother and proband WES sample for identification of variants related to a genetic disorder in the proband. The analysis method was adapted from the  Galaxy tutorial project <https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/exome-seq/tutorial.html>.

## Tools used

- BWA-MEM 2 (MAPPING)
- SAMBAMBA (BAM FILE PROCESSING)
- FREEBAYES ( JOINT VARIANT CALLING)
- BCFTOOLS (VCF FILE PROCESSING)
- SNPEFF (VARIANT ANNOTATION)
- VCFANNO (VARIANT ANNOTATION)

## Installation instructions

  1. Download Pixi: ` curl -fsSL https://pixi.sh/install.sh | sh `
  2. Add the following to ~/.bashrc: `eval "$(pixi completion --shell bash)"` to enable autocompletion. Restart the shell.
  3. Install the tools using Pixi: `pixi global install <tool_name>`
  4. Download the repo: ` git clone https://github.com/99ARITRA/Joint_genotyping.git `
  5. Navigate to the directory "Joint_genotyping" in terminal and run the script as shown below.

## Pipeline execution

 To execute the bash script type the following commands in a linux terminal or WSL terminal:
 `bash joint_genotyping.sh`
