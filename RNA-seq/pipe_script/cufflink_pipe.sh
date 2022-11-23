#!/bin/bash

# command line variable 
while getopts i:o:x:t: flag
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
    name=$(echo $file_name | grep -P '(?<=Tophat_result/)[a-zA-Z0-9]+' -o)
    echo $name # STAR
    cufflinks --num-threads ${thread} -o ${outdir}/CUFFLINKS/${name} -G ${indexdir} ${indir}/${name}/Aligned.sortedByCoord.out.bam
