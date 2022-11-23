#!/bin/bash

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



# find & unique, delemeter
file_list=$(find ${indir} -maxdepth 1 -type f -exec basename "{}" \; | cut -d'_' -f1 | sort -u)
echo $file_list

# run tophat2
for file_name in $file_list
do
    echo $file_name
    if [ ${type} = "fastq" ];then
        echo "Tophat2 Galore RUN!!!!!!!!"
        tophat2 --num-threads ${thread} -G ${indexdir}/gencode.v37.primary_assembly.annotation.gtf -o ${outdir}/Tophat_result/${file_name} ${indexdir}/hg38 ${indir}/${file_name}_R1_val_1.fq.gz ${indir}/${file_name}_R2_val_2.fq.gz
    else
        tophat2 --num-threads ${thread} -G ${indexdir}/gencode.v37.primary_assembly.annotation.gtf -o ${outdir}/Tophat_result/${file_name} ${indexdir}/hg38 ${indir}/${file_name}_R1.${type}.gz ${indir}/${file_name}_R2.${type}.gz
    fi
done
