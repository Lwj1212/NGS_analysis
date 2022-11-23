#!/bin/bash

# < QC + pipeline}
# 1. RNA_script/script/kallisto_salmon.sh -d [container_name] -r a -b [/home/wmbio/RNA_SEQ/temp/] -i [input_dir] -k [index_name] -t 20 -g fastq -c y
# < QC 이후, QC파일이 다른 폴더에 잇을 때>
# 2. RNA_script/script/kallisto_salmon.sh -d [container_name] -r a -b [/home/wmbio/RNA_SEQ/temp/] -i [qc_dir] -t 20 -g qc -c n

while getopts d:r:n:p:h:b:i:t:g:k:c: flag
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
        k) index=${OPTARG};;
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

# Salmon, Kallisto
echo "Trim Galore, SALMON, Kallisto"

docker run --rm -dit -v ${basedir}:/data --name ${containername}_trim_galore sempre813/trim_galore:0.6.6 /bin/bash
docker run --rm -dit -v ${basedir}:/data --name ${containername}_salmon combinelab/salmon:latest /bin/bash
docker run --rm -dit -v ${basedir}:/data --name ${containername}_kallisto zlskidmore/kallisto:latest /bin/bash

# dir making
docker exec -it ${containername}_trim_galore /bin/sh -c "mkdir -p /data/Paired_output_dir/QC; mkdir -p /data/Paired_output_dir/KALLISTO;mkdir -p /data/Paired_output_dir/SALMON;chmod -R 777 /data/Paired_output_dir;apt-get update && apt-get install -y pigz"

## RUN pipeline
# trim_galore CMD
if [ ${type} = "fastq" ];then
    echo "Trim Galore RUN!!!!!!!!"
    docker exec ${containername}_trim_galore /bin/sh -c "/data/NGS-Analysis/RNA-seq/pipe_script/trim_galore_pipe.sh -i /data/${indir} -o /data/Paired_output_dir"
fi

## RUN pipeline
if [ ${run} = "a" ];then

    if [ ${qc_indir} = "n" ];then
        # STAR CMD
        echo "Salmon RUN!!!!!!!!"
        docker exec ${containername}_salmon /bin/sh -c "/data/NGS-Analysis/RNA-seq/pipe_script/salmon_pipe.sh -i /data/${indir} -o /data/Paired_output_dir/SALMON -x /data/index/salmon_index/${index}.idx -t ${thread} -g ${type}"

        # kallisto CMD
        echo "Kallisto RUN!!!!!!!!"
        docker exec ${containername}_kallisto /bin/sh -c "/data/NGS-Analysis/RNA-seq/pipe_script/kallisto_pipe.sh -i /data/${indir} -o /data/Paired_output_dir/KALLISTO -x /data/index/kallisto_index/${index}.idx -t ${thread} -g ${type}"
    
    else
        # STAR CMD
        echo "Salmon RUN!!!!!!!!"
        docker exec ${containername}_salmon /bin/sh -c "/data/NGS-Analysis/RNA-seq/pipe_script/salmon_pipe.sh -i /data/Paired_output_dir/QC -o /data/Paired_output_dir/SALMON -x /data/index/salmon_index/${index}.idx -t ${thread} -g ${type}"

        # kallisto CMD
        echo "Kallisto RUN!!!!!!!!"
        docker exec ${containername}_kallisto /bin/sh -c "/data/NGS-Analysis/RNA-seq/pipe_script/kallisto_pipe.sh -i /data/Paired_output_dir/QC -o /data/Paired_output_dir/KALLISTO -x /data/index/kallisto_index/${index}.idx -t ${thread} -g ${type}"
    fi


elif [ ${run} = "k" ];then

    if [ ${qc_indir} = "n" ];then
        # kallisto CMD
        echo "Kallisto RUN!!!!!!!!"
        docker exec ${containername}_kallisto /bin/sh -c "/data/NGS-Analysis/RNA-seq/pipe_script/kallisto_pipe.sh -i /data/${indir} -o /data/Paired_output_dir/KALLISTO -x /data/index/kallisto_index/${index}.idx -t ${thread} -g ${type}"
    
    else
        # kallisto CMD
        echo "Kallisto RUN!!!!!!!!"
        docker exec ${containername}_kallisto /bin/sh -c "/data/NGS-Analysis/RNA-seq/pipe_script/kallisto_pipe.sh -i /data/Paired_output_dir/QC -o /data/Paired_output_dir/KALLISTO -x /data/index/kallisto_index/${index}.idx -t ${thread} -g ${type}"
    fi

elif [ ${run} = "s" ];then

    if [ ${qc_indir} = "n" ];then
        # STAR CMD
        echo "Salmon RUN!!!!!!!!"
        docker exec ${containername}_salmon /bin/sh -c "/data/NGS-Analysis/RNA-seq/pipe_script/salmon_pipe.sh -i /data/${indir} -o /data/Paired_output_dir/SALMON -x /data/index/salmon_index/${index}.idx -t ${thread} -g ${type}"
    else
        # STAR CMD
        echo "Salmon RUN!!!!!!!!"
        docker exec ${containername}_salmon /bin/sh -c "/data/NGS-Analysis/RNA-seq/pipe_script/salmon_pipe.sh -i /data/Paired_output_dir/QC -o /data/Paired_output_dir/SALMON -x /data/index/salmon_index/${index}.idx -t ${thread} -g ${type}"
    fi

else
    echo "No Run type!!!"
fi

docker stop ${containername}_kallisto ${containername}_trim_galore ${containername}_salmon
