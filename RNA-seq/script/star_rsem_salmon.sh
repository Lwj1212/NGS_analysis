#!/bin/bash

# ebi-browser example
# < QC + pipeline}
# 1. RNA_script/script/star_rsem_salmon.sh -r m -b /home/wmbio/RNA_SEQ/temp/ -i input_dir -t 30 -g fastq
# 2. RNA_script/script/star_rsem_salmon.sh -r a -b /home/wmbio/RNA_SEQ/temp/ -i input_dir -t 30 -g qc -s star_ensembl_RON_ADD -e rsem_ensembl_RON_ADD -c y

# < QC 이후, QC파일이 다른 폴더에 잇을 때, 반드시 -g qc, -c n 옵션을 사용할 것>
# <STAR, RSEM 모두 할 떄>
# RNA_script/script/star_rsem_salmon.sh -r a -b /home/wmbio/RNA_SEQ/temp/ -i QC_test -t 30 -g qc -s star_ensembl_RON_ADD -e rsem_ensembl_RON_ADD -c n
# < RSEM만 할 때, STAR의 폴더 위치를 지정>
# RNA_script/script/star_rsem_salmon.sh -r r -b /home/wmbio/RNA_SEQ/temp/ -i STAR -t 30 -g qc -s star_ensembl_RON_ADD -e rsem_ensembl_RON_ADD -c n

# index generate
# STAR --runThreadN 13 --runMode genomeGenerate --genomeDir /data/bak_star_rsem_index/star_index_RON_ADD --genomeFastaFiles /data/bak_star_rsem_index/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa  --sjdbGTFfile /data/bak_star_rsem_index/Homo_sapiens.GRCh38.104.chr.MST1R_del.sorted.gtf --sjdbOverhang 99
# rsem-prepare-reference --gtf /data/star_rsem_index/Homo_sapiens.GRCh38.104.chr.MST1R_del.sorted.gtf /data/star_rsem_index/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa rsem_index/rsem

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

# trimgalore
docker run --rm -dit -v ${basedir}:/data --name trim_galore sempre813/trim_galore:0.6.6 /bin/bash
docker run --rm -dit -v ${basedir}:/data --name star_rsem_pipe broadinstitute/gtex_rnaseq:V8 /bin/bash
docker run --rm -dit -v ${basedir}:/data --name salmon combinelab/salmon:latest /bin/bash

# dir making
docker exec -it star_rsem_pipe /bin/sh -c "mkdir -p /data/output_dir/QC; mkdir -p /data/output_dir/STAR;mkdir -p /data/output_dir/RSEM;mkdir -p /data/output_dir/SALMON;chmod -R 777 /data/output_dir"

## RUN pipeline
# trim_galore CMD
if [ ${type} = "fastq" ];then
    echo "Trim Galore RUN!!!!!!!!"
    docker exec trim_galore /bin/sh -c "/data/RNA_script/pipe_script/trim_galore_pipe.sh -i /data/${indir} -o /data/output_dir"
fi

## RUN pipeline
if [ ${run} = "a" ];then

    if [ ${qc_indir} = "n" ];then
        # STAR CMD
        echo "STAR RUN!!!!!!!!"
        docker exec star_rsem_pipe /bin/sh -c "/data/RNA_script/pipe_script/star_pipe.sh -i /data/${indir} -o /data/output_dir/STAR -x /data/index/star_index/${star_index} -t ${thread} -g ${type}"

        # RSEM CMD
        echo "RSEM RUN!!!!!!!!"
        docker exec star_rsem_pipe /bin/sh -c "/data/RNA_script/pipe_script/rsem_pipe.sh -i /data/output_dir/STAR -o /data/output_dir/RSEM -x /data/index/rsem_index/${rsem_index}/rsem -t ${thread}"

    else
        # STAR CMD
        echo "STAR RUN!!!!!!!!"
        docker exec star_rsem_pipe /bin/sh -c "/data/RNA_script/pipe_script/star_pipe.sh -i /data/output_dir/QC -o /data/output_dir/STAR -x /data/index/star_index/${star_index} -t ${thread} -g ${type}"

        # RSEM CMD
        echo "RSEM RUN!!!!!!!!"
        docker exec star_rsem_pipe /bin/sh -c "/data/RNA_script/pipe_script/rsem_pipe.sh -i /data/output_dir/STAR -o /data/output_dir/RSEM -x /data/index/rsem_index/${rsem_index}/rsem -t ${thread}"

    fi

elif [ ${run} = "s" ];then

    if [ ${qc_indir} = "n" ];then
        # STAR CMD
        echo "STAR RUN!!!!!!!!"
        docker exec star_rsem_pipe /bin/sh -c "/data/RNA_script/pipe_script/star_pipe.sh -i /data/${indir} -o /data/output_dir/STAR -x /data/index/star_index/${star_index} -t ${thread} -g ${type}"
    else
        # STAR CMD
        echo "STAR RUN!!!!!!!!"
        docker exec star_rsem_pipe /bin/sh -c "/data/RNA_script/pipe_script/star_pipe.sh -i /data/QC -o /data/output_dir/STAR -x /data/index/star_index/${star_index} -t ${thread} -g ${type}"
    fi

elif [ ${run} = "r" ];then

    if [ ${qc_indir} = "n" ];then
        # RSEM CMD
        echo "RSEM RUN!!!!!!!!"
        docker exec star_rsem_pipe /bin/sh -c "/data/RNA_script/pipe_script/rsem_pipe.sh -i /data/${indir} -o /data/output_dir/RSEM -x /data/index/rsem_index/${rsem_index}/rsem -t ${thread}"
    else
        # RSEM CMD
        echo "RSEM RUN!!!!!!!!"
        docker exec star_rsem_pipe /bin/sh -c "/data/RNA_script/pipe_script/rsem_pipe.sh -i /data/output_dir/STAR -o /data/output_dir/RSEM -x /data/index/rsem_index/${rsem_index}/rsem -t ${thread}"

    fi

elif [ ${run} = "m" ];then

    # salmon CMD
    echo "Salmon RUN!!!!!!!!"
    docker exec salmon /bin/sh -c "/data/RNA_script/pipe_script/salmon_pipe.sh -i /data/output_dir/QC -o /data/output_dir/SALMON -x /data/index/salmon_index/Homo_sapiens.GRCh38.cdna.all_add_RON_del -t ${thread} -g ${type}"

else
    echo "No Run type!!!"
fi

docker stop star_rsem_pipe trim_galore salmon