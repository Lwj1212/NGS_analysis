#!/bin/bash

while getopts i:o::x:t: flag
do
    case "${flag}" in
        i) indir=${OPTARG};;
        o) outdir=${OPTARG};;
        x) indexdir=${OPTARG};;
        t) thread=${OPTARG};;
    esac
done

# find & unique
for file_name in $indir/*
do
    echo $file_name
    name=$(echo $file_name | grep -P '(?<=STAR/)[a-zA-Z0-9]+' -o)
    echo $name
    mkdir -p ${outdir}/${name}
    rsem-calculate-expression --bam --no-bam-output -p ${thread} --paired-end ${indir}/${name}/Aligned.toTranscriptome.out.bam ${indexdir} ${outdir}/${name}/rsem
done
