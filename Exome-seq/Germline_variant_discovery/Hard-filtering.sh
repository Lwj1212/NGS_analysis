#!/bin/bash

START=$(date +%s)

while getopts v:r:b:o: flag
do
    case "${flag}" in
        v) mergevcf=${OPTARG};;
        r) ref=${OPTARG};;
        b) workhard=${OPTARG};;
        o) finalPath=${OPTARG};;
    esac
done

source /opt/conda/etc/profile.d/conda.sh
conda activate gatk4

# directory
mkdir -p ${finalPath}
mkdir -p ${workhard}

chmod 777 -R ${finalPath}

# SPLIT SNP
gatk --java-options '-Xmx25g' SelectVariants \
 -R ${ref} \
 -V ${mergevcf} \
 -select-type SNP \
 -O ${workhard}/Exome_Norm_HC_calls.snps.vcf

# SPLIT INDEL
gatk --java-options '-Xmx25g' SelectVariants \
 -R ${ref} \
 -V ${mergevcf} \
 -select-type INDEL \
 -O ${workhard}/Exome_Norm_HC_calls.indels.vcf

# FILTER SNP 
# GATK hard filter recommend
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
gatk --java-options '-Xmx25g' VariantFiltration \
 -R ${ref} \
 -V ${workhard}/Exome_Norm_HC_calls.snps.vcf \
 --missing-values-evaluate-as-failing true \
 --filter-expression "QD < 2.0" \
 --filter-name "QD2" \
 --filter-expression "QUAL < 30.0" \
 --filter-name "QUAL30" \
 --filter-expression "SOR > 3.0" \
 --filter-name "SOR3" \
 --filter-expression "FS > 60.0" \
 --filter-name "FS60" \
 --filter-expression "MQ < 40.0" \
 --filter-name "MQ40" \
 --filter-expression "MQRankSum < -12.5" \
 --filter-name "MQRankSum-12.5" \
 --filter-expression "ReadPosRankSum < -8.0" \
 --filter-name "ReadPosRankSum-8" \
 -O ${workhard}/HardFilter.snps.filtered.vcf

gatk --java-options '-Xmx25g' VariantFiltration \
 -R ${ref} \
 -V ${workhard}/Exome_Norm_HC_calls.indels.vcf \
 --missing-values-evaluate-as-failing true \
 --filter-expression "QD < 2.0" \
 --filter-name "QD2" \
 --filter-expression "QUAL < 30.0" \
 --filter-name "QUAL30" \
 --filter-expression "FS > 200.0" \
 --filter-name "FS200" \
 --filter-expression "ReadPosRankSum < -20.0" \
 --filter-name "ReadPosRankSum-20" \
 -O ${workhard}/HardFilter.indels.filtered.vcf

 # MERGE
gatk --java-options '-Xmx25g' MergeVcfs \
 -I ${workhard}/HardFilter.indels.filtered.vcf \
 -I ${workhard}/HardFilter.snps.filtered.vcf \
 -O ${finalPath}/HardFilter.filtered.vcf

# FILTER PASS ONLY
 gatk --java-options "-Xmx25g" SelectVariants \
  -R ${ref} \
  -V ${finalPath}/HardFilter.filtered.vcf \
  -O ${finalPath}/HardFilter.filtered.PASS.vcf \
  --exclude-filtered

# Time stemp
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "HardFilter $DIFF seconds"