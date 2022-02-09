#!/bin/bash

SAMPLE=$1

for SNVS in $(ls *$SAMPLE*/results/variants/somatic.snvs.vcf.gz)
do
  DIR_NAME=$(echo $SNVS | awk -F"/" '{print $1}')

  NEW_SNVS=$DIR_NAME".somatic.snvs.vcf.gz"
  NEW_SNVS=${NEW_SNVS/.somatic.somatic/.somatic}

  echo "$SNVS -> $NEW_SNVS"
  cp $SNVS $NEW_SNVS

  INDELS=${SNVS/snvs.vcf.gz/indels.vcf.gz}
  NEW_INDELS=$DIR_NAME".somatic.indels.vcf.gz"
  NEW_INDELS=${NEW_INDELS/.somatic.somatic/.somatic}

  echo "$INDELS -> $NEW_INDELS"
  cp $INDELS $NEW_INDELS

  tar cvzpf $DIR_NAME".strelka2s.tgz" $DIR_NAME
done
