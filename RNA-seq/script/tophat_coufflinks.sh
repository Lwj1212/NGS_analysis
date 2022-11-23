#!/bin/bash

# tophat, cufflink
# RUN ex) RNA_script/script/[.sh] -r a -n jinoo -p Wmlswkdia1! -h 192.168.0.90 -b /home/wmbio/temp -i input_dir -x index -t 30 -g fq
# command line variable 
while getopts r:n:p:h:b:i:t:g: flag
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
docker stop trim_galore tophat2 cufflink
echo "Trim Galore, Tophat2, Cufflink, Docker Container RUN!!!"
# trimgalore
# docker run --rm -dit -v ${basedir}:/data --name trim_galore dukegcb/trim-galore:latest /bin/bash
docker run --rm -dit -v ${basedir}:/data --name trim_galore sempre813/trim_galore:0.6.6 /bin/bash
# tophat container
docker run --rm -dit -v ${basedir}:/data --name tophat2 genomicpariscentre/tophat2:latest /bin/bash
# cufflink container
docker run --rm -dit -v ${basedir}:/data --name cufflink fomightez/rnaseqcufflinks:latest /bin/bash

# dir making
docker exec -it tophat2 /bin/sh -c "mkdir -p /data/output_dir/QC; mkdir -p /data/output_dir/Tophat_result;mkdir -p /data/output_dir/cufflinks_result;chmod -R 777 /data/output_dir;apt-get update && apt-get install -y pigz"

## RUN pipeline
# trim_galore CMD
if [ ${type} = "fastq" ];then
    echo "Trim Galore RUN!!!!!!!!"
    docker exec trim_galore /bin/sh -c "/data/RNA_script/pipe_script/trim_galore_pipe.sh -i /data/${indir} -o /data/output_dir"
else
    docker exec -it trim_galore /bin/sh -c "cp /data/${indir}/* /data/output_dir/QC/"
fi

if [ ${run} = "a" ];then

    # tophat CMD
    echo "Tophat2 RUN!!!!!!!!"
    docker exec tophat2 /bin/sh -c "/data/RNA_script/pipe_script/tophat_pipe.sh -i /data/output_dir/QC -o /data/output_dir -x /data/index/genome_index -t ${thread} -g ${type}"

    # cufflink CMD
    echo "Cufflinks RUN!!!!!!!!"
    docker exec cufflink /bin/sh -c "/data/RNA_script/pipe_script/cufflink_pipe.sh -i /data/output_dir/Tophat_result -o /data/output_dir -x /data/index/genome_index -t ${thread}"
    # done
    echo "DONE!!!"
else
    echo "No Run type!!!"
fi

docker stop trim_galore tophat2 cufflink