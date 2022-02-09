# CRISPR-toolbox
Version: 2019-03-27
Author: Taejoon Kwon (tkwon@unist.ac.kr)

## Introduction

This is the document of how to design crRNA to attack multiple InDels in cancer cells. 

## Step 0. Mapping
Currently, I am mapping paired FASTQ files in the following procedure.

1. Rename the files as "<Project Name>_<Cell Type>_<Unique sample ID>_R[12].raw.fastq" 
     (i.e. SNU201903_Gliblastoma+e_GBL-67-tumor_R1.raw.fastq). 
     "+e" is used for an exome-sequencing.
     Use "-tumor" or "-normal" at the end of 'Unique sample ID'.

2. Trim the reads with trimmomatic (/scripts/run-trimmomatic.sh). 
     Install trimmomatic using conda ('conda -c bioconda install trimmomatic).

3. Map the reads with BWA MEM, and make 'sort-index-rmdup-indexed' BAM file 
     using samtools (/scripts/run-bwa_p.sh).

## Step 1. Variant Calling

1. Make a config directory to run Strelka2 (/scripts/conf-strelka2-somatic.sh).
2. Run a Strelka2 (/scripts/run-strelka2.sh).
3. Rename VCF file and compress the file (/scripts/extract-strelka2-vcf.sh).

## Step 2. Variant Filtering

1. Filter InDel from VCF with 3-8 length range 
  * /vcf/filter-vcf-InDel.py
2. Filter by Read Depth.
  * /vcf/filter-vcf-DP.py for Strelka2 somatic variants
  * /vcf/filter-vcf-AD.py for Strelka2 germline variants
3. Filter by variant scores (optional)
  * /vcf/filter-vcf-strelka2SomaticEVC.py 
4. Filter by homozygosity (optional)
  * /vcf/filter-vcf-Homo.py

## Step 3. crRNA Design

1. Make var_flank_seq from vcf and reference fasta.
  * /crRNA/make-var_flank_seq.py <filtered vcf>

2. Make crRNA fasta from var_flank_seq.
  * /crRNA/var_flank_seq-to-crRNA <.var_flank_seq>

3. Make alternative short fragments for all variants.
  * /crRNA/make-var_flank_seq.py <Strelka2 raw vcf without filtering>
  * /crRNA/var_flank_seq-to-alt_fasta.py

4. Run exonerate to check the off-targets
  * exonerate <FASTA for designed crRNAs> <ref genome>
  * exonerate <FASTA for designed crRNAs> <alt seq from Step 3.3 above>

5. Select crRNAs without off target candidates
  * /crRNA/select-crRNA.py <crRNA FASTA> <exonerate out to ref> <exonerate out to alt>

## Step 4. Report 

1. Copy /report/make-region_tview_script.py to the directory, and run it with
  designed crRNA FASTA file. It will generate bash commands to run samtools tview.
  * /report/make-region_tview_script.py <crRNA FASTA> > run.sh
  * Check the location of BAM files and the referene FASTA file.

2. Run the script produced by the above step. It will generate tview files for each 
  crRNA target site.
  * bash run.sh

3. Make the summary report by running /report/make-somatic-tview-report.txt. 
