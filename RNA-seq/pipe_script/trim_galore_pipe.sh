#!/bin/bash

while getopts i:o::x:t:g: flag
do
    case "${flag}" in
        i) indir=${OPTARG};;
        o) outdir=${OPTARG};;
    esac
done

# find & unique, delemeter
file_list=$(find ${indir} -maxdepth 1 -type f -exec basename "{}" \; | cut -d'_' -f1 | sort -u)
echo $file_list

for file_name in $file_list
do
    echo $file_name  #K1Aligned_Access_R1.fastq.gz
    trim_galore -j 4 --no_report_file --q 30 --paired --clip_R1 10 --clip_R2 10 ${indir}/${file_name}_R1.fastq.gz ${indir}/${file_name}_R2.fastq.gz -o ${outdir}/QC 
done
