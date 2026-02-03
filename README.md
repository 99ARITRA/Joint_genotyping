# Joint genotyping
This is a project on joint genotyping using NGS tools on a father, mother and proband sample for identification of variants related to Osteopertrosis. The analysis method was adapted from a Galaxy tutorial project on the disease.
## Tools
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
  - GENOTYPEGVCFSS
- VCF PROCESSING
  - BCFTOOLS NORM
  - BCFTOOLS FILLTAGS
- VARIANT FILTRATION
  - GATK VARIANTFILTRATION
 ## Pipeline execution
 To execute the bash script type the following commands in a linux terminal or WSL terminal:
 `bash germline_VC.sh`
