#!/bin/bash
NUM_THREADS=8

DB="/home/taejoon/pub/db.bwa/HUMAN_ens90_dna_rm"
DBNAME=$(basename $DB)

for FQ1 in $(ls ../fastq/*_1P)
do
  FQ2=${FQ1/_1P/_2P}

  SAM=$(basename $FQ1)
  SAM=${SAM/_1P/}"."$DBNAME".bwa_mem.sam"
  BAM=${SAM/.sam/}".bam"
  SORTED_BAM=${SAM/.sam/}".sorted.bam"

  if [ -e $SORTED_BAM ]; then
    echo "Skip $FQ1"
    continue
  fi

  if [ ! -e $SAM ]; then
    echo "Make $SAM"
    bwa mem -t $NUM_THREADS $DB $FQ1 $FQ2 > $SAM
  fi
    

  if [ ! -e $BAM ]; then
    echo "Make $BAM"
    samtools view -h -@ $NUM_THREADS -o $BAM $SAM
  fi

  if [ -e $BAM ] && [ ! -e $SORTED_BAM ]; then
    echo "Make $SORTED_BAM"
    samtools sort -@ $NUM_THREADS -o $SORTED_BAM $BAM
    samtools index -@ $NUM_THREADS $SORTED_BAM
  fi

  RMDUP_BAM=${SAM/.sam/}".rmdup.bam"
  if [ -e $SORTED_BAM ] && [ ! -e $RMDUP_BAM ]; then
    echo "Make $RMDUP_BAM"
    samtools rmdup -s --output-fmt BAM $SORTED_BAM $RMDUP_BAM
  fi
  
  BAI=$RMDUP_BAM".bai"
  if [ -e $RMDUP_BAM ] && [ ! -e $BAI ]; then
    echo "Make $BAI "
    samtools index -@ $NUM_THREADS $RMDUP_BAM
  fi
done
