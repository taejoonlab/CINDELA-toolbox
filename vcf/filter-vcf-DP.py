#!/usr/bin/env python3
# Filter VCF records based on minimum "Read Depth (DP)"
import gzip
import sys

filename_vcf = sys.argv[1]
min_DP_cutoff = int(sys.argv[2])

filename_out = filename_vcf.replace('.vcf', '')

f_vcf = open(filename_vcf, 'r')
if filename_vcf.endswith('.gz'):
    f_vcf = gzip.open(filename_vcf, 'rt')
    filename_out = filename_vcf.replace('.vcf.gz', '')

filename_out += ('+DP%d.vcf' % min_DP_cutoff)
f_out = open(filename_out, 'w')
for line in f_vcf:
    if line.startswith('#'):
        f_out.write(line.strip() + "\n")
        continue

    tokens = line.strip().split("\t")
    tmp_ref = tokens[3]
    tmp_alt = tokens[4]
    tmp_fmt_tokens = tokens[8].split(':')
    idx_DP = tmp_fmt_tokens.index('DP')

    min_DP = -1
    for tmp in tokens[9:]:
        tmp_DP = sum(int(x) for x in tmp.split(':')[idx_DP].split(','))
        if min_DP < 0:
            min_DP = tmp_DP
        elif tmp_DP < min_DP:
            min_DP = tmp_DP

    if min_DP > min_DP_cutoff:
        f_out.write(line.strip() + "\n")
f_out.close()
