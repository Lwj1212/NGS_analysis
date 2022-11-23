#!/bin/bash

source /opt/conda/etc/profile.d/conda.sh

# <samtool env>
conda create -n samtools
conda activate samtools
conda install -c bioconda samtools bwa

# <gatk4 env>
conda create -n gatk4
conda activate gatk4
conda install -c bioconda gatk4=4.1.8.1 bcftools

# <vep env>
conda create -n vep
conda activate vep
conda install -c bioconda ensembl-vep=105.0 samtools
