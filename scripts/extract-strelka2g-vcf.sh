#!/bin/bash

SAMPLE=$1

for VCF in $(ls *$SAMPLE*/results/variants/variants.vcf.gz)
do
  DIR_NAME=$(echo $VCF | awk -F"/" '{print $1}')
  NEW=$DIR_NAME".germline.variants.vcf.gz"
  NEW=${NEW/.germline.germline/.germline}
  echo "$VCF -> $NEW"

  cp $VCF $NEW
  tar cvzpf $DIR_NAME".strelka2g.tgz" $DIR_NAME
done
