#!/usr/bin/env python3
import gzip
import sys

filename_vcf = sys.argv[1]
filename_filter_vcf = sys.argv[2]
filter_name = sys.argv[3]

filename_out = filename_vcf.replace('.vcf', '')

filter_list = dict()
f_filter = open(filename_filter_vcf, 'r')
if filename_filter_vcf.endswith('.gz'):
    f_filter = gzip.open(filename_filter_vcf, 'rt')
    filename_out = filename_vcf.replace('.vcf.gz', '')

for line in f_filter:
    if line.startswith('#'):
        continue
    tokens = line.strip().split("\t")
    var_id = '%s:%s:%s-%s' % (tokens[0], tokens[1], tokens[3], tokens[4])
    var_name = tokens[2]
    filter_list[var_id] = var_name
f_filter.close()

f_vcf = open(filename_vcf, 'r')
if filename_vcf.endswith('.gz'):
    f_vcf = gzip.open(filename_vcf, 'rt')

count_total = 0
count_filtered = 0
f_out = open('%s.+No%s.vcf' % (filename_out, filter_name), 'w')
f_drop = open('%s.+No%s.drop.vcf' % (filename_out, filter_name), 'w')
for line in f_vcf:
    if line.startswith('#'):
        f_out.write("%s\n" % line.strip())
        f_drop.write("%s\n" % line.strip())
        continue

    tokens = line.strip().split("\t")
    tmp_ref = tokens[3]
    tmp_alt = tokens[4]
    tmp_var_id = '%s:%s:%s-%s' % (tokens[0], tokens[1], tokens[3], tokens[4])
    count_total += 1

    if tokens[0].startswith('Mm'):
        continue

    if tmp_var_id not in filter_list:
        count_filtered += 1
        f_out.write("%s\n" % line.strip())
    else:
        f_drop.write("%s\n" % line.strip())
f_out.close()
f_drop.close()

sys.stderr.write('Total Vars: %d\n' % count_total)
sys.stderr.write('Filtered Vars: %d\n' % count_filtered)
