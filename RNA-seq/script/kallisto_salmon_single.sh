#!/bin/bash

# ebi-browser example
# 1. RNA_script/script/[.sh] -r a -b /home/wmbio/ -i PUBLIC_DATA -t 15 -g fastq

# index generate
# STAR --runThreadN 13 --runMode genomeGenerate --genomeDir /data/bak_star_rsem_index/star_index_RON_ADD --genomeFastaFiles /data/bak_star_rsem_index/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa  --sjdbGTFfile /data/bak_star_rsem_index/Homo_sapiens.GRCh38.104.chr.MST1R_del.sorted.gtf --sjdbOverhang 99
# rsem-prepare-reference --gtf /data/star_rsem_index/Homo_sapiens.GRCh38.104.chr.MST1R_del.sorted.gtf /data/star_rsem_index/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa rsem_index/rsem

while getopts d:r:n:p:h:b:i:t:g:s:e:c: flag
do
    case "${flag}" in
        d) containername=${OPTARG};;
        r) run=${OPTARG};;
        n) nasid=${OPTARG};;
        p) naspw=${OPTARG};;
        h) nas=${OPTARG};;
        b) basedir=${OPTARG};;
        i) indir=${OPTARG};;
        t) thread=${OPTARG};;
        g) type=${OPTARG};;
        s) star_index=${OPTARG};;
        e) rsem_index=${OPTARG};;
        c) qc_indir=${OPTARG};;
    esac
done

# index file download
if [ ! -d ${basedir}/index ];then
    echo "Index dir not exist!"
    ncftpget -u ${nasid} -p ${naspw} ${nas} ${basedir} /Data/RNA_SEQ/work/index.tar.gz
    tar -xzvf index.tar.gz -C ${basedir}
else
    echo "Index dir exist!"
fi

## Docker run
docker stop kallisto trim_galore salmon && docker rm kallisto trim_galore salmon

# Salmon, Kallisto
echo "Trim Galore, SALMON, Kallisto"

docker run --rm -dit -v ${basedir}:/data --name ${containername}_trim_galore sempre813/trim_galore:0.6.6 /bin/bash
docker run --rm -dit -v ${basedir}:/data --name ${containername}_salmon combinelab/salmon:latest /bin/bash
docker run --rm -dit -v ${basedir}:/data --name ${containername}_kallisto zlskidmore/kallisto:latest /bin/bash

# dir making
docker exec -it ${containername}_trim_galore /bin/sh -c "mkdir -p /data/Single_output_dir/QC; mkdir -p /data/Single_output_dir/KALLISTO;mkdir -p /data/Single_output_dir/SALMON;chmod -R 777 /data/Single_output_dir;apt-get update && apt-get install -y pigz"

## RUN pipeline
# trim_galore CMD
if [ ${type} = "fastq" ];then
    echo "Trim Galore RUN!!!!!!!!"
    docker exec ${containername}_trim_galore /bin/sh -c "/data/RNA_script/pipe_script/trim_galore_pipe_single.sh -i /data/${indir} -o /data/Single_output_dir"
fi

## RUN pipeline
if [ ${run} = "a" ];then
    if [ ${qc_indir} = "n" ];then
        echo "Kallisto RUN!!!!!!!!"
        docker exec ${containername}_kallisto /bin/sh -c "/data/RNA_script/pipe_script/kallisto_pipe_single.sh -i /data/${indir} -o /data/Single_output_dir/KALLISTO -x /data/index/kallisto_index/Homo_sapiens.GRCh38.cdna.all_add_RON_del.idx -t ${thread} -g ${type}"

        echo "Salmon RUN!!!!!!!!"
        docker exec ${containername}_salmon /bin/sh -c "/data/RNA_script/pipe_script/salmon_pipe_single.sh -i /data/${indir} -o /data/Single_output_dir/SALMON -x /data/index/salmon_index/Homo_sapiens.GRCh38.cdna.all_add_RON_del -t ${thread} -g ${type}"

    else
        echo "Kallisto RUN!!!!!!!!"
        docker exec ${containername}_kallisto /bin/sh -c "/data/RNA_script/pipe_script/kallisto_pipe_single.sh -i /data/Single_output_dir/QC -o /data/Single_output_dir/KALLISTO -x /data/index/kallisto_index/Homo_sapiens.GRCh38.cdna.all_add_RON_del.idx -t ${thread} -g ${type}"

        echo "Salmon RUN!!!!!!!!"
        docker exec ${containername}_salmon /bin/sh -c "/data/RNA_script/pipe_script/salmon_pipe_single.sh -i /data/Single_output_dir/QC -o /data/Single_output_dir/SALMON -x /data/index/salmon_index/Homo_sapiens.GRCh38.cdna.all_add_RON_del -t ${thread} -g ${type}"

    fi

elif [ ${run} = "k" ];then

    if [ ${qc_indir} = "n" ];then
        # kallisto CMD
        echo "Kallisto RUN!!!!!!!!"
        docker exec ${containername}_kallisto /bin/sh -c "/data/RNA_script/pipe_script/kallisto_pipe_single.sh -i /data/${indir} -o /data/Single_output_dir/KALLISTO -x /data/index/kallisto_index/Homo_sapiens.GRCh38.cdna.all_add_RON_del.idx -t ${thread} -g ${type}"

    else
        # kallisto CMD
        echo "Kallisto RUN!!!!!!!!"
        docker exec ${containername}_kallisto /bin/sh -c "/data/RNA_script/pipe_script/kallisto_pipe_single.sh -i /data/Single_output_dir/QC -o /data/Single_output_dir/KALLISTO -x /data/index/kallisto_index/Homo_sapiens.GRCh38.cdna.all_add_RON_del.idx -t ${thread} -g ${type}"
    fi

elif [ ${run} = "s" ];then
    if [ ${qc_indir} = "n" ];then
        echo "Salmon RUN!!!!!!!!"
        docker exec ${containername}_salmon /bin/sh -c "/data/RNA_script/pipe_script/salmon_pipe_single.sh -i /data/${indir} -o /data/Single_output_dir/SALMON -x /data/index/salmon_index/Homo_sapiens.GRCh38.cdna.all_add_RON_del -t ${thread} -g ${type}"
    else
        echo "Salmon RUN!!!!!!!!"
        docker exec ${containername}_salmon /bin/sh -c "/data/RNA_script/pipe_script/salmon_pipe_single.sh -i /data/Single_output_dir/QC -o /data/Single_output_dir/SALMON -x /data/index/salmon_index/Homo_sapiens.GRCh38.cdna.all_add_RON_del -t ${thread} -g ${type}"
    fi
else
    echo "No Run type!!!"
fi

docker stop ${containername}_kallisto ${containername}_trim_galore ${containername}_salmon
