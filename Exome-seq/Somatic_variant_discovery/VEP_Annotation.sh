#!/bin/bash

START=$(date +%s)

while getopts c:r:o: flag
do
    case "${flag}" in
        c) vcfFolder=${OPTARG};;
        r) ref=${OPTARG};;
        o) finalPath=${OPTARG};;
    esac
done

source /opt/conda/etc/profile.d/conda.sh
conda activate vep

mkdir -p ${finalPath}

# VEP and VCF2MAF
for mapFile in ${vcfFolder}/*.mt2_filtered_sort_pass.vcf
do

    filename=$(basename $mapFile .mt2_filtered_sort_pass.vcf)
    output=${filename}.vep.maf

    perl /root/vcf2maf-1.6.21/vcf2maf.pl \
        --input-vcf ${mapFile} \
		--output-maf ${finalPath}/${output} \
        --tumor-id ${filename} \
		--vep-path /opt/vep/src/ensembl-vep \
		--vep-data /root/work/vep_resource/ \
        --vep-overwrite \
		--ref-fasta ${ref} \
		--ncbi-build GRCh38
done

# Time stemp
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "Filter_Annotation Using Funcotation $DIFF seconds"