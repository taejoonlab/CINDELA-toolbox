#!/usr/bin/env python3
import sys
import gzip

filename_1 = sys.argv[1]
filename_2 = sys.argv[2]

sample_1 = sys.argv[3]
sample_2 = sys.argv[4]


def read_vcf(tmp_filename):
    rv = dict()
    f = open(tmp_filename, 'r')
    if tmp_filename.endswith('.gz'):
        f = gzip.open(tmp_filename, 'rt')

    for line in f:
        if line.startswith('#'):
            continue
        tokens = line.strip().split("\t")
        tmp_seq_id = tokens[0]
        tmp_pos = tokens[1]
        tmp_ref = tokens[3]
        tmp_alt = tokens[4]
        tmp_tag = '%s|%s|%s-%s' % (tmp_seq_id, tmp_pos, tmp_ref, tmp_alt)
        rv[tmp_tag] = 1
    f.close()
    return rv


rv_1 = read_vcf(filename_1)
rv_2 = read_vcf(filename_2)
set_1 = set(rv_1.keys())
set_2 = set(rv_2.keys())

print("%d = overlap" % len(set_1.intersection(set_2)))
print("%d = %s only" % (len(set_1 - set_2), filename_1))
print("%d = %s only" % (len(set_2 - set_1), filename_2))

output_name = '%s-%s' % (sample_1, sample_2)
f_overlap = open('%s.overlap.vcf' % output_name, 'w')

for tmp in list(set_1.intersection(set_2)):
    f_overlap.write('%s\n' % tmp)
f_overlap.close()

f_only_1 = open('%s.%s_only.vcf' % (output_name, sample_1), 'w')
for tmp in list(set_1 - set_2):
    f_only_1.write('%s\n' % tmp)
f_only_1.close()

f_only_2 = open('%s.%s_only.vcf' % (output_name, sample_2), 'w')
for tmp in list(set_2 - set_1):
    f_only_2.write('%s\n' % tmp)
f_only_2.close()
