#!/bin/bash

# ebi-browser example
# 1. RNA_script/script/[.sh]  -r m -b /home/wmbio/ -i PUBLIC_DATA -t 15 -g fastq
# 2. RNA_script/script/[.sh] -r a -b /home/wmbio/ -i PUBLIC-2-5 -t 15 -g qc -s star_ensembl_RON_ADD -e rsem_ensembl_RON_ADD

# index generate
# STAR --runThreadN 13 --runMode genomeGenerate --genomeDir /data/bak_star_rsem_index/star_index_RON_ADD --genomeFastaFiles /data/bak_star_rsem_index/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa  --sjdbGTFfile /data/bak_star_rsem_index/Homo_sapiens.GRCh38.104.chr.MST1R_del.sorted.gtf --sjdbOverhang 99
# rsem-prepare-reference --gtf /data/star_rsem_index/Homo_sapiens.GRCh38.104.chr.MST1R_del.sorted.gtf /data/star_rsem_index/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa rsem_index/rsem

while getopts r:n:p:h:b:i:t:g:s:e: flag
do
    case "${flag}" in
        r) run=${OPTARG};;
        n) nasid=${OPTARG};;
        p) naspw=${OPTARG};;
        h) nas=${OPTARG};;
        b) basedir=${OPTARG};;
        i) indir=${OPTARG};;
        t) thread=${OPTARG};;
        g) type=${OPTARG};;
        s) star_index=${OPTARG};;
        e) cufflink_index=${OPTARG};;
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
docker stop star_rsem_pipe trim_galore && docker rm star_rsem_pipe trim_galore

# trimgalore
docker run --rm -dit -v ${basedir}:/data --name trim_galore sempre813/trim_galore:0.6.6 /bin/bash
docker run --rm -dit -v ${basedir}:/data --name star_rsem_pipe broadinstitute/gtex_rnaseq:V8 /bin/bash
docker run --rm -dit -v ${basedir}:/data --name cufflink fomightez/rnaseqcufflinks:latest /bin/bash

# dir making
docker exec -it star_rsem_pipe /bin/sh -c "mkdir -p /data/output_dir/QC; mkdir -p /data/output_dir/STAR;mkdir -p /data/output_dir/CUFFLINKS;chmod -R 777 /data/output_dir;apt-get update && apt-get install -y pigz"

## RUN pipeline
# trim_galore CMD
if [ ${type} = "fastq" ];then
    echo "Trim Galore RUN!!!!!!!!"
    docker exec trim_galore /bin/sh -c "/data/RNA_script/pipe_script/trim_galore_pipe.sh -i /data/${indir} -o /data/output_dir"
fi

## RUN pipeline
if [ ${run} = "a" ];then
    # STAR CMD
    echo "STAR RUN!!!!!!!!"
    docker exec star_rsem_pipe /bin/sh -c "/data/RNA_script/pipe_script/star_pipe.sh -i /data/output_dir/QC -o /data/output_dir/STAR -x /data/index/star_index/${star_index} -t ${thread} -g ${type}"

    # cufflink CMD
    echo "Cufflinks RUN!!!!!!!!"
    docker exec cufflink /bin/sh -c "/data/RNA_script/pipe_script/cufflink_pipe.sh -i /data/output_dir/STAR -o /data/output_dir -x /data/index/genome_index/${cufflink_index} -t ${thread}"
 
elif [ ${run} = "s" ];then
    # STAR CMD
    echo "STAR RUN!!!!!!!!"
    docker exec star_rsem_pipe /bin/sh -c "/data/RNA_script/pipe_script/star_pipe.sh -i /data/output_dir/QC -o /data/output_dir/STAR -x /data/index/star_index/${star_index} -t ${thread} -g ${type}"

elif [ ${run} = "c" ];then
    # cufflink CMD
    echo "Cufflinks RUN!!!!!!!!"
    docker exec cufflink /bin/sh -c "/data/RNA_script/pipe_script/cufflink_pipe.sh -i /data/output_dir/STAR -o /data/output_dir -x /data/index/genome_index/${cufflink_index} -t ${thread}"

else
    echo "No Run type!!!"
fi

docker stop star_rsem_pipe trim_galore cufflink