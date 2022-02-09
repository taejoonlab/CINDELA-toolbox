#!/usr/bin/env python3
import gzip
import sys

filename_vcf = sys.argv[1]

filename_out = filename_vcf.replace('.vcf', '')
f_vcf = open(filename_vcf, 'r')
if filename_vcf.endswith('.gz'):
    f_vcf = gzip.open(filename_vcf, 'rt')
    filename_out = filename_vcf.replace('.vcf.gz', '')

filename_out += '.+PASS.vcf'
f_out = open(filename_out, 'w')

count_total = 0
count_pass = 0
count_chrUn = 0

for line in f_vcf:
    if line.startswith('#'):
        f_out.write(line.strip()+"\n")
        continue

    tokens = line.strip().split("\t")
    tmp_chr = tokens[0]

    tmp_ref = tokens[3]
    tmp_alt = tokens[4]

    tmp_filter = tokens[6]

    count_total += 1
    # Skip without PASS call
    if tmp_filter.find('PASS') < 0:
        continue

    if tmp_chr.startswith('chrUn') or tmp_chr.find('_random') >= 0:
        count_chrUn += 1
        continue

    count_pass += 1
    f_out.write(line.strip()+"\n")

f_vcf.close()
f_out.close()

sys.stderr.write('Total variants: %d\n' % count_total)
sys.stderr.write('Variants on chrUn: %d\n' % count_chrUn)
sys.stderr.write('Passed variants: %d\n' % count_pass)
