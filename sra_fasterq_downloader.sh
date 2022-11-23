#!/bin/bash

while getopts b:f:i: flag
do
    case "${flag}" in
        b) basedir=${OPTARG};;
        f) foldername=${OPTARG};;
        i) input=${OPTARG};;
    esac
done
 
while IFS= read -r line || [ -n "$line" ]; do

    echo $line
    docker run -t --rm -v ${basedir}/${foldername}:/output:rw -w /output ncbi/sra-tools fasterq-dump -e 15 -p --split-files ${line}

done < $input

# # file rename
sudo chmod 777 ${basedir}/${foldername}
python ${basedir}/script/ena_file_rename.py ${basedir}/${foldername}
gzip ${basedir}/${foldername}/*