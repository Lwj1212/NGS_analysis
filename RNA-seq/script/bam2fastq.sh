#!/bin/bash

# In samtools17 container

dir=/raw_data/RNASeq_Access/bam
outdir=/raw_data/RNASeq_Access/bamtofastq
workdir=/raw_data/RNASeq_Access/work
name=$1 # mother dir

for bam_name in $dir/*sortedByCoord.out.bam
do
    echo $bam_name
    file_name=$(echo $bam_name | grep -P '(?<=bam/)[a-zA-Z0-9]+' -o)
    echo $file_name

    ## samtools
    samtools view --threads 15 -u -f 1 -F 12 ${bam_name} > ${workdir}/${file_name}_map_map.bam
    # R1 unmapped, R2 mapped
    samtools view --threads 15 -u -f 4 -F 264 ${bam_name} > ${workdir}/${file_name}_unmap_map.bam
    # R1 mapped, R2 unmapped
    samtools view --threads 15 -u -f 8 -F 260 ${bam_name} > ${workdir}/${file_name}_map_unmap.bam
    # R1 & R2 unmapped
    samtools view --threads 15 -u -f 12 -F 256 ${bam_name} > ${workdir}/${file_name}_unmap_unmap.bam

    # merge
    samtools merge --threads 15 -u ${workdir}/${file_name}_unmapped.bam ${workdir}/${file_name}_unmap_map.bam ${workdir}/${file_name}_map_unmap.bam ${workdir}/${file_name}_unmap_unmap.bam

    # sort
    samtools sort --threads 15 -n ${workdir}/${file_name}_map_map.bam -o ${workdir}/${file_name}_mapped.sort.bam
    samtools sort --threads 15 -n ${workdir}/${file_name}_unmapped.bam -o ${workdir}/${file_name}_unmapped.sort.bam

    ## bedtools
    echo "bamtofastq"
    bamToFastq -i ${workdir}/${file_name}_mapped.sort.bam -fq ${workdir}/${file_name}_mapped.1.fastq -fq2 ${workdir}/${file_name}_mapped.2.fastq 
    bamToFastq -i ${workdir}/${file_name}_unmapped.sort.bam -fq ${workdir}/${file_name}_unmapped.1.fastq -fq2 ${workdir}/${file_name}_unmapped.2.fastq 

    # betools merge
    cat ${workdir}/${file_name}_mapped.1.fastq ${workdir}/${file_name}_unmapped.1.fastq > ${outdir}/${file_name}_Access_R1.fastq
    cat ${workdir}/${file_name}_mapped.2.fastq ${workdir}/${file_name}_unmapped.2.fastq > ${outdir}/${file_name}_Access_R2.fastq

done