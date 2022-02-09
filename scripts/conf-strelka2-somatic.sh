#!/bin/bash
REF="$HOME/pub/ucsc/hg38/hg38_chr.fa"

# Setup the Strelka2
# https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/quickStart.md
# 
# Current version: 2.9.10
# https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2

DIR_STRELKA="$HOME/src/strelka/2.9.10/"

BAM_TUMOR="name of tumor BAM file"
BAM_NORMAL="name of normal BAM file"

OUT=${BAM_TUMOR/.bam/}".somatic"
OUT=$(basename $OUT)
echo $OUT

# add '--exome' for Exome sequencing
$DIR_STRELKA/bin/configureStrelkaSomaticWorkflow.py --referenceFasta=$REF --runDir=$OUT \
       --normalBam=$BAM_NORMAL --tumorBam=$BAM_TUMOR

