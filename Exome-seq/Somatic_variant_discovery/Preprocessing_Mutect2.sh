#!/bin/bash
#!/bin/bash

START=$(date +%s)

while getopts v:r:g:i:b: flag
do
    case "${flag}" in
        v) bamFolder=${OPTARG};;
        r) ref=${OPTARG};;
        g) gnomad=${OPTARG};;
        i) interval=${OPTARG};;
        b) workvariant=${OPTARG};;
    esac
done

source /opt/conda/etc/profile.d/conda.sh
conda activate gatk4

# directory
mkdir -p ${workvariant}

# Tumor-only somatic variant call pipeline
for mapFile in ${bamFolder}/*_final.bam
do  
    filename=$(basename ${mapFile} _bwa_final.bam)

    # 1. Generate OXOG metrics:
    gatk --java-options "-Xmx25G -Xmx25G" CollectSequencingArtifactMetrics \
        -I ${mapFile} \
        -O ${workvariant}/${filename} --FILE_EXTENSION .txt \
        -R ${ref}

    # 2. Generate pileup summaries on tumor sample
    # interval
    for i in `seq -f %04g 0 14`
    do  
        outfile=${filename}_targeted_sequencing_${i}.table

        gatk GetPileupSummaries \
        -I ${mapFile} \
        -V ${gnomad} \
        -L ${interval}/${i}-scattered.interval_list \
        -R ${ref} \
        -O ${workvariant}/${outfile} &
    done
    wait

    # 3. Merge pileup summaries merge
    head -2 ${workvariant}/${filename}_targeted_sequencing_0000.table > ${workvariant}/${filename}_targeted_sequencing.table
    tail -n +3 -q ${workvariant}/${filename}_targeted_sequencing_00* >> ${workvariant}/${filename}_targeted_sequencing.table

    ## 4. Calculate contamination on tumor sample
    infile=${filename}_targeted_sequencing.table
    outfile=${filename}_targeted_sequencing.contamination.table

    gatk CalculateContamination \
        -I ${workvariant}/${infile} \
        -O ${workvariant}/${outfile}

    ## 5. Find tumor sample name from BAM
    gatk GetSampleName \
        -I ${mapFile} \
        -O ${workvariant}/${filename}.targeted_sequencing.sample_name

done


# Time stemp
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "Preprocessing SomaticVariant Calling $DIFF seconds"