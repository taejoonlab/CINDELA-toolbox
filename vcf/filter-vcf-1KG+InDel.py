#!/usr/bin/env python3
import os
import gzip
import sys

filename_vcf = sys.argv[1]

var_list = dict()

dirname_1kg = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           '..', 'data')
filename_1kg = os.path.join(dirname_1kg, '1000g_pon.hg38.+InDel2-8.vcf.gz')
f_1kg = gzip.open(filename_1kg, 'rt')
sys.stderr.write('Read %s\n' % filename_1kg)
for line in f_1kg:
    if line.startswith('#'):
        continue
    tokens = line.strip().split("\t")
    for tmp_alt in tokens[4].split(','):
        var_id = '%s:%s:%s-%s' % (tokens[0], tokens[1], tokens[3], tmp_alt)
        var_name = tokens[2]
        var_list[var_id] = var_name
f_1kg.close()


f_vcf = open(filename_vcf, 'r')
if filename_vcf.endswith('.gz'):
    f_vcf = gzip.open(filename_vcf, 'rt')

count_total = 0
count_filtered = 0
filename_out = filename_vcf.replace('.vcf', '+No1KG.vcf')
if filename_vcf.endswith('.gz'):
    filename_out = filename_vcf.replace('.vcf.gz', '+No1KG.vcf')

f_out = open(filename_out, 'w')
for line in f_vcf:
    if line.startswith('#'):
        f_out.write("%s\n" % line.strip())
        continue

    tokens = line.strip().split("\t")
    tmp_ref = tokens[3]
    tmp_alt = tokens[4]
    tmp_var_id = '%s:%s:%s-%s' % (tokens[0], tokens[1], tokens[3], tokens[4])
    count_total += 1

    if tokens[0].startswith('Mm'):
        continue

    if tmp_var_id not in var_list:
        count_filtered += 1
        f_out.write("%s\n" % line.strip())
f_out.close()

sys.stderr.write('Total Vars: %d\n' % count_total)
sys.stderr.write('Filtered Vars: %d\n' % count_filtered)
