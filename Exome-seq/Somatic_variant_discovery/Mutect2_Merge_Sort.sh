#!/bin/bash

START=$(date +%s)

while getopts v:r:d:g:p:i:b:o: flag
do
    case "${flag}" in
        v) bamFolder=${OPTARG};;
        r) ref=${OPTARG};;
        d) refDict=${OPTARG};;
        g) gnomad=${OPTARG};;
        p) pon=${OPTARG};;
        i) interval=${OPTARG};;
        b) workvariant=${OPTARG};;
        o) workVCF=${OPTARG};;
    esac
done

source /opt/conda/etc/profile.d/conda.sh
conda activate gatk4

mkdir -p ${workVCF}

## Run MuTect2 using only tumor sample on chromosome level (25 commands with different intervals)    
for mapFile in ${bamFolder}/*_final.bam
do
    for i in `seq -f %04g 0 14`
    do  

    filename=$(basename $mapFile _bwa_final.bam)
    output=${filename}.mt2_${i}.vcf
    flr2_output=${filename}.flr2_${i}.tar.gz

    # file list
    if [ ${i} = 1 ]
    then
        echo ${workVCF}/${output} > ${workVCF}/${filename}_m2_vcf_file.list

    else
        echo ${workVCF}/${output} >> ${workVCF}/${filename}_m2_vcf_file.list
    fi

    gatk --java-options "-Xmx33G -XX:ConcGCThreads=1" Mutect2 \
        -R ${ref} \
        -L ${interval}/${i}-scattered.interval_list \
        -I ${mapFile} \
        --native-pair-hmm-threads 2 \
        -tumor ${workvariant}/${filename}.targeted_sequencing.sample_name \
        --af-of-alleles-not-in-resource 2.5e-06 \
        --germline-resource ${gnomad} \
        --f1r2-tar-gz ${workVCF}/${flr2_output} \
        -pon ${pon} \
        --tmp-dir /root \
        -O ${workVCF}/${output} &
    done
    wait

    # merge scattered phenotype vcf files & filter
    filename=$(basename $mapFile _bwa_final.bam)
    combine=${filename}.mt2_merged.vcf
    sort=${filename}.mt2_merged_sort.vcf

    # Merge
    gatk --java-options "-Xmx20G" GatherVcfs \
            -R ${ref} \
            -I ${workVCF}/${filename}_m2_vcf_file.list \
            -O ${workVCF}/${combine}

    # Sort
    gatk --java-options "-Xmx20G" SortVcf \
            --SEQUENCE_DICTIONARY ${refDict} \
            --CREATE_INDEX true \
            -I ${workVCF}/${combine} \
            -O ${workVCF}/${sort}


    # Stats Merge
    statslist=$(for f in ${workVCF}/${filename}*.mt2_*.vcf.stats; do echo -n "--stats $f " ;done)
    gatk MergeMutectStats \
        ${statslist} \
        -O ${workVCF}/${sort}.stats

    # Flr2 Merge
    flr2list=$(for f in ${workVCF}/${filename}*.flr2_*.tar.gz; do echo -n "-I $f " ;done)
    gatk LearnReadOrientationModel \
        ${flr2list} \
        -O ${workVCF}/${filename}_read-orientation-model.tar.gz


done

# Time stemp
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "Mutect2 Merge $DIFF seconds"