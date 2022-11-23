#!/bin/bash

while getopts r:d:t: flag
do
    case "${flag}" in
        r) ref=${OPTARG};;
        d) dir=${OPTARG};;
        t) threads=${OPTARG};;
    esac
done

source /opt/conda/etc/profile.d/conda.sh
conda activate gatk4

# ScatterIntervalsByNs [Picard type]
gatk ScatterIntervalsByNs -R ${ref} -OT ACGT -O temp.interval_list

# Splitintervals --> 가용한 Core에 맞춰 짜를 것.
gatk SplitIntervals -R ${ref} -L temp.interval_list --scatter-count ${threads} -O ${dir}