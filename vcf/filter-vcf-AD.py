#!/usr/bin/env python3
# Filter VCF records based on minimum "Allelic Depth (AD)"
import gzip
import sys

filename_vcf = sys.argv[1]
min_AD_cutoff = int(sys.argv[2])
# min_AD_cutoff = 15

filename_out = filename_vcf.replace('.vcf', '')

f_vcf = open(filename_vcf, 'r')
if filename_vcf.endswith('.gz'):
    f_vcf = gzip.open(filename_vcf, 'rt')
    filename_out = filename_vcf.replace('.vcf.gz', '')

count_var = 0
filename_out += ('.+AD%d.vcf' % min_AD_cutoff)
f_out = open(filename_out, 'w')
for line in f_vcf:
    if line.startswith('#'):
        f_out.write(line.strip() + "\n")
        continue

    tokens = line.strip().split("\t")
    tmp_ref = tokens[3]
    tmp_alt = tokens[4]
    tmp_fmt_tokens = tokens[8].split(':')
    idx_AD = tmp_fmt_tokens.index('AD')

    min_AD = -1
    for tmp in tokens[9:]:
        tmp_AD = sum([int(x) for x in tmp.split(':')[idx_AD].split(',')])
        if min_AD < 0:
            min_AD = tmp_AD
        elif tmp_AD < min_AD:
            min_AD = tmp_AD

    if min_AD > min_AD_cutoff:
        f_out.write(line.strip() + "\n")
        count_var += 1
f_out.close()
sys.stderr.write('Count Var: %d\n' % count_var)
