#!/usr/bin/env python3
import gzip
import sys

filename_vcf = sys.argv[1]

filename_out = filename_vcf.replace('.vcf', '')
f_vcf = open(filename_vcf, 'r')
if filename_vcf.endswith('.gz'):
    f_vcf = gzip.open(filename_vcf, 'rt')
    filename_out = filename_vcf.replace('.vcf.gz', '')

filename_out += '.+Homo.vcf'
f_out = open(filename_out, 'w')

count_total = 0
count_filtered = 0
for line in f_vcf:
    if line.startswith('#'):
        f_out.write(line.strip()+"\n")
        continue

    tokens = line.strip().split("\t")
    tmp_ref = tokens[3]
    tmp_alt = tokens[4]
    count_total += 1

    tmp_gtype_alt = tokens[9].split(':')[0].count('1')

    # only homo
    if tmp_gtype_alt == 2:
        f_out.write(line.strip()+"\n")
        count_filtered += 1
f_out.close()

sys.stderr.write('Total vars: %d\n' % count_total)
sys.stderr.write('Filtered vars: %d\n' % count_filtered)
