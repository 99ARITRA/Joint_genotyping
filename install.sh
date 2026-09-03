#!/usr/bin/bash
set -e
# Update the image
apt update -y && apt upgrade -y
apt install wget -y
apt install bsdextrautils -y
# Install pixi
wget -qO- https://pixi.sh/install.sh | bash
export PATH="/root/.pixi/bin:$PATH"
echo 'export PATH="/root/.pixi/bin:$PATH"' >> /root/.bashrc
source /root/.bashrc
# Add autocompletion for pixi
echo -e "eval \"\$(pixi completion --shell bash)\"" >> /root/.bashrc
# Enable bioconda
pixi config set default-channels '["conda-forge", "bioconda"]'
# Install required packages
pixi global install --environment ngs bwa-mem2 sambamba freebayes bcftools snpeff snpsift vcfanno gatk4
pixi global install --environment py python pandas
apt install tabix -y