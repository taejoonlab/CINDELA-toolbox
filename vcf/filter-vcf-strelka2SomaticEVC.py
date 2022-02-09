#!/usr/bin/env python3
# Filter VCF records based on "Somatic Empirical Variant Score(SomaticEVS)"
# reported by Strelka2
import gzip
import sys

filename_vcf = sys.argv[1]
min_SEVS_cutoff = int(sys.argv[2])

filename_out = filename_vcf.replace('.vcf', '')

f_vcf = open(filename_vcf, 'r')
if filename_vcf.endswith('.gz'):
    f_vcf = gzip.open(filename_vcf, 'rt')
    filename_out = filename_vcf.replace('.vcf.gz', '')


def parse_record(tmp_str):
    rv = dict()
    for tmp in tmp_str.split(';')[1:]:
        if tmp.find('=') < 0:
            continue
        (tmp_k, tmp_v) = tmp.split('=')
        rv[tmp_k] = tmp_v
    return rv


count_var = 0
filename_out += ('+SomaticEVS%d.vcf' % min_SEVS_cutoff)
f_out = open(filename_out, 'w')
for line in f_vcf:
    if line.startswith('#'):
        f_out.write(line.strip() + "\n")
        continue

    tokens = line.strip().split("\t")
    tmp_ref = tokens[3]
    tmp_alt = tokens[4]
    tmp_rec = parse_record(tokens[7])

    if float(tmp_rec['SomaticEVS']) > min_SEVS_cutoff:
        f_out.write(line.strip() + "\n")
        count_var += 1
f_out.close()

sys.stderr.write('Filtered: %d\n' % count_var)
