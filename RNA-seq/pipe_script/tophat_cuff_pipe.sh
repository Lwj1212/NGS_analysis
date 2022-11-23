#!/bin/bash

# tophat, cufflink
# RUN ex) tophat_cuff_pipe.sh -n NASid -p NASpw -h NASip -b /home/wmbio/temp -i input_dir -x index -t 10 -g fq
# command line variable 
while getopts n:p:h:b:i:x:t:g: flag
do
    case "${flag}" in
        n) nasid=${OPTARG};;
        p) naspw=${OPTARG};;
        h) nas=${OPTARG};;
        b) basedir=${OPTARG};;
        i) indir=${OPTARG};;
        x) indexdir=${OPTARG};;
        t) thread=${OPTARG};;
        g) type=${OPTARG};;
    esac
done

# # RNA-Seq script
# git clone https://github.com/Jin0331/RNA_script.git ${basedir}/RNA_script


# index file download
if [ ! -d ${basedir}/index ];then
    echo "No Index dir"
    ncftpget -u ${nasid} -p ${naspw} ${nas} ${basedir} /Data/RNA_SEQ/work/index.tar.gz
    tar -xzvf index.tar.gz -C ${basedir}
else
    echo "Index dir exist!"
fi


## Docker run
# trimgalore
docker run --rm -dit -v ${basedir}:/data --name trim_galore dukegcb/trim-galore:latest /bin/bash

# tophat container
docker run --rm -dit -v ${basedir}:/data --name tophat2 genomicpariscentre/tophat2:latest /bin/bash

# cufflink container
docker run --rm -dit -v ${basedir}:/data --name cufflink fomightez/rnaseqcufflinks:latest /bin/bash

# dir making
docker exec -it tophat2 /bin/sh -c "mkdir -p /data/output_dir/QC; mkdir -p /data/output_dir/Tophat_result;mkdir -p /data/output_dir/cufflinks_result;chmod -R 777 /data/output_dir"

# trim_galore CMD
echo "Trim Galore RUN!!!!!!!!"
docker exec trim_galore /bin/sh -c "/data/RNA_script/trim_galore_pipe.sh -i /data/${indir} -o /data/output_dir"

# tophat CMD
echo "Tophat2 RUN!!!!!!!!"
docker exec tophat2 /bin/sh -c "/data/RNA_script/tophat_pipe.sh -i /data/output_dir/QC -o /data/output_dir -x /data/${indexdir}/genome_index -t ${thread} -g ${type}"

# cufflink CMD
echo "CUFFLINKS RUN!!!!!!!!"
docker exec cufflink /bin/sh -c "/data/RNA_script/cufflink_pipe.sh -i /data/output_dir/Tophat_result -o /data/output_dir -x /data/${indexdir}/genome_index -t ${thread}"

# done
docker stop trim_galore tophat2 cufflink