#!/bin/bash

CR_FA=$1

if [ ! -e $CR_FA ]; then
  echo "CR_FA is not available."
  exit
fi

if [ ! -e REF.fa ]; then
  echo "REF.fa is not available."
  exit
else
  OUT=${CR_FA/.fa/}".REF.vulgar"

  REF_INFO=$(readlink REF.fa)
  echo "#REF: " $REF_INFO
  echo "#REF: " $REF_INFO > $OUT

  echo "Run exonerate ... "$CR_FA" to REF.fa"
  exonerate --bestn 5 $CR_FA REF.fa | grep ^vulgar >> $OUT
fi

if [ ! -e SOMATIC.fa ]; then
  echo "SOMATIC.fa is not available"
else
  OUT=${CR_FA/.fa/}".SOMATIC.vulgar"

  ALT_INFO=$(readlink SOMATIC.fa)
  echo "#SOMATIC: " $ALT_INFO
  echo "#SOMATIC: " $ALT_INFO > $OUT

  echo "Run exonerate ... "$CR_FA" to SOMATIC.fa"
  exonerate --bestn 50 $CR_FA SOMATIC.fa | grep ^vulgar >> $OUT
fi

if [ ! -e GERMLINE.fa ]; then
  echo "GERMLINE.fa is not available"
else
  OUT=${CR_FA/.fa/}".GERMLINE.vulgar"

  ALT_INFO=$(readlink GERMLINE.fa)
  echo "#GERMLINE: " $ALT_INFO
  echo "#GERMLINE: " $ALT_INFO > $OUT

  echo "Run exonerate ... "$CR_FA" to GERMLINE.fa"
  exonerate --bestn 5 $CR_FA GERMLINE.fa | grep ^vulgar >> $OUT
fi
