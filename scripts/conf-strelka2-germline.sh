#!/bin/bash

# Make a symbolic link to reference genome

REF="$HOME/pub/ucsc/hg38/hg38_chr.fa"

# Setup the Strelka2
# https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/quickStart.md
# 
# Current version: 2.9.10
# https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2

DIR_STRELKA="$HOME/src/strelka/2.9.10/"

for BAM in $(ls ../bwa/*.bam)
do
  OUT=${BAM/.bam/}".germline"
  OUT=$(basename $OUT)
  echo $OUT

  # add '--exome' for Exome sequencing
  # add '--rna' for RNA-seq
  $DIR_STRELKA/bin/configureStrelkaGermlineWorkflow.py --referenceFasta=$REF --runDir=$OUT --bam=$BAM
done
