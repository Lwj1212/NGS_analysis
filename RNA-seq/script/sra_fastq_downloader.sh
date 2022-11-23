#!/bin/bash

# argument
while getopts b:i:c:f: flag
do
    case "${flag}" in
        b) basedir=${OPTARG};;
        i) filename=${OPTARG};;
        c) containername=${OPTARG};;
        f) foldername=${OPTARG};;
    esac
done
 
docker pull pegi3s/sratoolkit 
while read line; do
# reading each line
    echo $line
    sudo docker run -dit -v ${basedir}/:/data --name ${containername}_${line} pegi3s/sratoolkit /bin/bash
    
    sudo docker exec ${containername}_${line} /bin/bash -c "mkdir /data/${foldername}"
    sudo docker exec ${containername}_${line} fastq-dump --gzip -F --split-files ${line} --outdir /data/${foldername}
    sudo docker exec ${containername}_${line} /bin/bash -c "chmod -R 777 /data/${foldername}"
    sudo docker exec ${containername}_${line} /bin/bash -c "mv /root/ncbi_error_report.xml /data/"
    sudo docker stop ${containername}_${line}
    sudo docker rm ${containername}_${line}

    mv ${basedir}/ncbi_error_report.xml ${basedir}/${line}_ncbi_error_report.xml
done < $filename

# file rename
python ${basedir}/RNA_script/script/ena_file_rename.py ${basedir}/${foldername}
