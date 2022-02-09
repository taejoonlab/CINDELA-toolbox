#!/bin/bash

## Assume that you install trimmomatic with conda
## Directory for conda

DIR_CONDA=$HOME/miniconda3
NUM_THREADS=3

FA_ADAPTER=$( find $DIR_CONDA | grep -m1 TruSeq3-PE.fa )
echo "Adapter: "$FA_ADAPTER
cp $FA_ADAPTER .

pigz -p 4 -d *fastq.gz

for FQ1 in $(ls *R1.*fastq)
do
  FQ2=${FQ1/_R1/_R2}
  OUT=${FQ1/_R1.raw.fastq}"_trim"
  echo $FQ1 $FQ2 $OUT

  trimmomatic PE -threads $NUM_THREADS -validatePairs \
   -summary $OUT".summary" $FQ1 $FQ2 -baseout $OUT \
   ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
done

