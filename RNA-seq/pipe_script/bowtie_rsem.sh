#!/bin/bash

# broadinstitute/gtex_rnaseq:V8
# apt update && apt install bowtie2  ---> /usr/bin에 bowtie2 설치됨.

while getopts i:o::x:t:g: flag
do
    case "${flag}" in
        i) indir=${OPTARG};;
        o) outdir=${OPTARG};;
        x) indexdir=${OPTARG};;
        t) thread=${OPTARG};;
        g) type=${OPTARG};;
    esac
done

file_list=$(find ${indir} -maxdepth 1 -type f -exec basename "{}" \; | cut -d'_' -f1 | sort -u)
echo $file_list

# run tophat2
for file_name in $file_list
do
    echo $file_name  #K1Aligned_Access_R1.fastq.gz
    mkdir -p ${outdir}/${file_name}
    rsem-calculate-expression -p ${thread} --paired-end --bowtie2 --bowtie2-path /usr/bin --estimate-rspd --no-bam-output --append-names ${indir}/${file_name}_R1.${type}.gz ${indir}/${file_name}_R2.${type}.gz ${indexdir} ${outdir}/${file_name}/${file_name}
done


