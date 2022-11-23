#!/bin/bash

START=$(date +%s)

while getopts v:r:b:h:n:k:d:m:a:o: flag
do
    case "${flag}" in
        v) mergevcf=${OPTARG};;
        r) ref=${OPTARG};;
        b) workvqsr=${OPTARG};;
        h) hapmap=${OPTARG};;
        n) omni=${OPTARG};;
        k) kg1000=${OPTARG};;
        d) dbsnp=${OPTARG};;
        m) mills=${OPTARG};;
        a) axiom=${OPTARG};;
        o) finalPath=${OPTARG};;
    esac
done

source /opt/conda/etc/profile.d/conda.sh
conda activate gatk4

# directory
mkdir -p ${workvqsr}
mkdir -p ${finalPath}

chmod -R 777 ${workvqsr}
chmod -R 777 ${finalPath}

# Variant quality score recalibration - SNPs
gatk --java-options "-Xms4G -Xmx4G -XX:ParallelGCThreads=5" VariantRecalibrator \
  -mode SNP \
  -R ${ref} \
  -V ${mergevcf} \
  -O ${workvqsr}/recalibrate_SNP.recal \
  --trust-all-polymorphic \
  -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
  --max-gaussians 6 \
  -tranche 100.0 -tranche 99.95 -tranche 99.9 \
  -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
  -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 \
  -tranche 92.0 -tranche 91.0 -tranche 90.0 \
  --tranches-file ${workvqsr}/recalibrate_SNP.tranches \
  --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap} \
  --resource:omni,known=false,training=true,truth=false,prior=12.0 ${omni} \
  --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${kg1000} \
  --resource:dbsnp,known=true,training=false,truth=false,prior=7.0 ${dbsnp} \
  --rscript-file ${workvqsr}/output_SNP1.plots.R

# Variant quality score recalibration - INDELs
gatk --java-options "-Xms4G -Xmx4G -XX:ParallelGCThreads=5" VariantRecalibrator \
  -mode INDEL \
  -R ${ref} \
  -V ${mergevcf} \
  -O ${workvqsr}/recalibrate_INDEL.recal \
  --trust-all-polymorphic \
  -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
  --max-gaussians 4 \
  -tranche 100.0 -tranche 99.95 -tranche 99.9 \
  -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
  -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 \
  -tranche 92.0 -tranche 91.0 -tranche 90.0 \
  --tranches-file ${workvqsr}/recalibrate_INDEL.tranches \
  --resource:mills,known=false,training=true,truth=true,prior=12.0 ${mills} \
  # --resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${axiom} \
  --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp} \
  --rscript-file ${workvqsr}/output_indel1.plots.R

  # Apply recalibration to SNPs
gatk --java-options "-Xms4G -Xmx4G -XX:ParallelGCThreads=5" ApplyVQSR \
  -V ${mergevcf} \
  -O ${workvqsr}/recalibrated_snps_VQSR.vcf \
  --recal-file ${workvqsr}/recalibrate_SNP.recal \
  --tranches-file ${workvqsr}/recalibrate_SNP.tranches \
  -truth-sensitivity-filter-level 99.7 \
  --create-output-variant-index true \
  -mode SNP

  # Apply recalibration to Indels
gatk --java-options "-Xms4G -Xmx4G -XX:ParallelGCThreads=5" ApplyVQSR \
  -V ${workvqsr}/recalibrated_snps_VQSR.vcf \
  -O ${finalPath}/recalibrated_VQSR.vcf \
  --recal-file ${workvqsr}/recalibrate_INDEL.recal \
  --tranches-file ${workvqsr}/recalibrate_INDEL.tranches \
  -truth-sensitivity-filter-level 99.7 \
  --create-output-variant-index true \
  -mode INDEL

gatk --java-options "-Xmx25g" SelectVariants \
  -R ${ref} \
  -V ${finalPath}/recalibrated_VQSR.vcf \
  -O ${finalPath}/recalibrated_VQSR.filtered.PASS.vcf \
  --exclude-filtered

# Time stemp
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "VQSR $DIFF seconds"