#!/usr/bin/env python3
import os
import gzip
import sys

filename_vcf = sys.argv[1]

gnomAD_list = dict()

dirname_gnomAD = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                              '..', 'data')

for tmp_filename in os.listdir(dirname_gnomAD):
    if not tmp_filename.endswith('vcf.gz'):
        continue
    if tmp_filename.find('gnomad') < 0:
        continue

    filename_gnomAD = os.path.join(dirname_gnomAD, tmp_filename)
    f_gnomAD = gzip.open(filename_gnomAD, 'rt')
    sys.stderr.write('Read %s\n' % filename_gnomAD)
    for line in f_gnomAD:
        tokens = line.strip().split("\t")
        for tmp_alt in tokens[4].split(','):
            var_id = '%s:%s:%s-%s' % (tokens[0], tokens[1], tokens[3], tmp_alt)
            var_name = tokens[2]
            gnomAD_list[var_id] = var_name
    f_gnomAD.close()

f_vcf = open(filename_vcf, 'r')
if filename_vcf.endswith('.gz'):
    f_vcf = gzip.open(filename_vcf, 'rt')

count_total = 0
count_filtered = 0
filename_out = filename_vcf.replace('.vcf', '+NoGnomAD.vcf')
if filename_vcf.endswith('.gz'):
    filename_out = filename_vcf.replace('.vcf.gz', '+NoGnomAD.vcf')

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

    if tmp_var_id not in gnomAD_list:
        count_filtered += 1
        f_out.write("%s\n" % line.strip())
f_out.close()

sys.stderr.write('Total Vars: %d\n' % count_total)
sys.stderr.write('Filtered Vars: %d\n' % count_filtered)
