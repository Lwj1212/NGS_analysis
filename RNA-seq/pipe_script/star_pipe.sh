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

# run STAR
for file_name in $file_list
do
    echo $file_name
    mkdir -p ${outdir}/${file_name}
    if [ ${type} = "fastq" ];then
        STAR --genomeDir ${indexdir} --readFilesIn ${indir}/${file_name}_R1_val_1.fq.gz ${indir}/${file_name}_R2_val_2.fq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic --runThreadN ${thread} --outFileNamePrefix ${outdir}/${file_name}/
    elif [ ${type} = "fq" ];then
        STAR --genomeDir ${indexdir} --readFilesIn ${indir}/${file_name}_R1.fq.gz ${indir}/${file_name}_R2.fq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic --runThreadN ${thread} --outFileNamePrefix ${outdir}/${file_name}/
    else
        STAR --genomeDir ${indexdir} --readFilesIn ${indir}/${file_name}_R1_val_1.fq.gz ${indir}/${file_name}_R2_val_2.fq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic --runThreadN ${thread} --outFileNamePrefix ${outdir}/${file_name}/
    fi
done
