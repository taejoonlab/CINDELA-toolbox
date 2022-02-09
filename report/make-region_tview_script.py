#!/usr/bin/env python3
import os
import sys

filename_ref = 'REF.fa'
filename_bam_list = 'BAM_LIST'
    
f_bam_list = open(filename_bam_list, 'r')
bam_list = f_bam_list.readlines()
f_bam_list.close()

filename_fa = sys.argv[1]

f_cr = open(filename_fa, 'r')
target_name = filename_fa.split('.')[1]

cr_idx = 1
for line in f_cr:
    if not line.startswith('>'):
        continue

    tokens = line.strip().lstrip('>').split("|")
    tmp_cr_id = tokens[0]
    tmp_var_id = tokens[1]
    tmp_chr = tmp_var_id.split(':')[0]
    tmp_pos = int(tmp_var_id.split(':')[1]) - 30
    tmp_region = "%s:%d" % (tmp_chr, tmp_pos)
    
    for filename_bam in bam_list:
        filename_bam = filename_bam.strip()
        tmp_bam_name = os.path.basename(filename_bam).split('.')[0]
        filename_out = '%s.%s.%s' % (tmp_cr_id, tmp_var_id.replace(':','_'), tmp_bam_name)
        print("echo %s" %filename_out)
        #print("samtools tview ../%s -d text -p %s > %s.dna_raw" % (filename_bam, tmp_region, filename_out))
        print("samtools tview %s -d text -p %s %s> %s.tview" % (filename_bam, tmp_region, filename_ref, filename_out))
f_cr.close()
