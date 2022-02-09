#!/bin/bash

MSI_LIST="$HOME/pub/ucsc/hg38/hg38_ref.MSIsensor-pro.list"

BAM_N="<BAM for normal sample>"
BAM_T="<BAM for tumor sample>"

OUT=$(basename $BAM_T)
OUT=${OUT/.bam/}".MSIsensor"

# Install MSIsensor-pro via conda
# $ conda install -c bioconda msisensor-pro

msisensor-pro msi -d $MSI_LIST -n $BAM_N -t $BAM_T -o $OUT
