#!/bin/bash

# ebi-browser example
# 1. RNA_script/script/wmbio_pipe_2_single.sh -r m -b /home/wmbio/ -i PUBLIC_DATA -t 15 -g fastq
# RUN ex) RNA_script/script/wmbio_pipe_2_single.sh -r a -b /home/wmbio/ -i PUBLIC-2-5 -t 15 -g qc -s star_ensembl_RON_ADD -e rsem_ensembl_RON_ADD

while getopts r:n:p:h:b:i:t:g:s:e:c: flag
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
        e) rsem_index=${OPTARG};;
        c) qc_indir=${OPTARG};;
    esac
done

# # RNA-Seq script
# git clone https://github.com/Jin0331/RNA_script.git ${basedir}/RNA_script


# index file download
if [ ! -d ${basedir}/index ];then
    echo "Index dir not exist!"
    ncftpget -u ${nasid} -p ${naspw} ${nas} ${basedir} /Data/RNA_SEQ/work/index.tar.gz
    tar -xzvf index.tar.gz -C ${basedir}
else
    echo "Index dir exist!"
fi

## Docker run
docker stop star_rsem_pipe trim_galore salmon && docker rm star_rsem_pipe trim_galore salmon

# Broad Institue STAR + RSEM PIPELINE DOCKER IMAGE
echo "Trim Galore, STAR + RSEM Container RUN!!!"

# trimgalore
docker run --rm -dit -v ${basedir}:/data --name trim_galore sempre813/trim_galore:0.6.6 /bin/bash
docker run --rm -dit -v ${basedir}:/data --name star_rsem_pipe broadinstitute/gtex_rnaseq:V8 /bin/bash
docker run --rm -dit -v ${basedir}:/data --name salmon combinelab/salmon:latest /bin/bash

# dir making
docker exec -it star_rsem_pipe /bin/sh -c "mkdir -p /data/output_dir/QC; mkdir -p /data/output_dir/STAR;mkdir -p /data/output_dir/RSEM;mkdir -p /data/output_dir/SALMON;chmod -R 777 /data/output_dir"

## RUN pipeline
if [ ${run} = "a" ];then

    if [ ${qc_indir} = "n" ];then
        # STAR CMD
        echo "STAR RUN!!!!!!!!"
        docker exec star_rsem_pipe /bin/sh -c "/data/RNA_script/pipe_script/star_pipe_single.sh -i /data/${indir} -o /data/output_dir/STAR -x /data/index/star_index/${star_index} -t ${thread} -g ${type}"

        # RSEM CMD
        echo "RSEM RUN!!!!!!!!"
        docker exec star_rsem_pipe /bin/sh -c "/data/RNA_script/pipe_script/rsem_pipe_single.sh -i /data/output_dir/STAR -o /data/output_dir/RSEM -x /data/index/rsem_index/${rsem_index}/rsem -t ${thread}"

    else
        # STAR CMD
        echo "STAR RUN!!!!!!!!"
        docker exec star_rsem_pipe /bin/sh -c "/data/RNA_script/pipe_script/star_pipe_single.sh -i /data/output_dir/QC -o /data/output_dir/STAR -x /data/index/star_index/${star_index} -t ${thread} -g ${type}"

        # RSEM CMD
        echo "RSEM RUN!!!!!!!!"
        docker exec star_rsem_pipe /bin/sh -c "/data/RNA_script/pipe_script/rsem_pipe_single.sh -i /data/output_dir/STAR -o /data/output_dir/RSEM -x /data/index/rsem_index/${rsem_index}/rsem -t ${thread}"

    fi

elif [ ${run} = "s" ];then

    if [ ${qc_indir} = "n" ];then
        # STAR CMD
        echo "STAR RUN!!!!!!!!"
        docker exec star_rsem_pipe /bin/sh -c "/data/RNA_script/pipe_script/star_pipe_single.sh -i /data/${indir} -o /data/output_dir/STAR -x /data/index/star_index/${star_index} -t ${thread} -g ${type}"
    else
        # STAR CMD
        echo "STAR RUN!!!!!!!!"
        docker exec star_rsem_pipe /bin/sh -c "/data/RNA_script/pipe_script/star_pipe_single.sh -i /data/QC -o /data/output_dir/STAR -x /data/index/star_index/${star_index} -t ${thread} -g ${type}"
    fi

elif [ ${run} = "r" ];then

    if [ ${qc_indir} = "n" ];then
        # RSEM CMD
        echo "RSEM RUN!!!!!!!!"
        docker exec star_rsem_pipe /bin/sh -c "/data/RNA_script/pipe_script/rsem_pipe_single.sh -i /data/${indir} -o /data/output_dir/RSEM -x /data/index/rsem_index/${rsem_index}/rsem -t ${thread}"
    else
        # RSEM CMD
        echo "RSEM RUN!!!!!!!!"
        docker exec star_rsem_pipe /bin/sh -c "/data/RNA_script/pipe_script/rsem_pipe_single.sh -i /data/output_dir/STAR -o /data/output_dir/RSEM -x /data/index/rsem_index/${rsem_index}/rsem -t ${thread}"

    fi

elif [ ${run} = "m" ];then

    # salmon CMD
    echo "Salmon RUN!!!!!!!!"
    docker exec salmon /bin/sh -c "/data/RNA_script/pipe_script/salmon_pipe_single.sh -i /data/output_dir/QC -o /data/output_dir/SALMON -x /data/index/salmon_index/Homo_sapiens.GRCh38.cdna.all_add_RON_del -t ${thread} -g ${type}"

else
    echo "No Run type!!!"
fi

docker stop star_rsem_pipe trim_galore salmon