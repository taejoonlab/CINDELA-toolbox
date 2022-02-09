#!/usr/bin/env python3
import gzip
import sys

filename_vcf = sys.argv[1]

filename_masked_fa = '/home/taejoon/pub/hg38/hg38.ucsc/hg38_masked.fa.gz'

seq_list = dict()
f_masked = open(filename_masked_fa, 'r')
if filename_masked_fa.endswith('.gz'):
    f_masked = gzip.open(filename_masked_fa, 'rt')
for line in f_masked:
    if line.startswith('>'):
        tmp_h = line.strip().lstrip('>')
        seq_list[tmp_h] = []
    else:
        seq_list[tmp_h].append(line.strip())
f_masked.close()


seq_dict = dict()
for tmp_h in seq_list.keys():
    seq_dict[tmp_h] = ''.join(seq_list[tmp_h])

f_vcf = open(filename_vcf, 'r')
filename_out = filename_vcf.replace('.vcf', '')
f_vcf = open(filename_vcf, 'r')
if filename_vcf.endswith('.gz'):
    f_vcf = gzip.open(filename_vcf, 'rt')
    filename_out = filename_vcf.replace('.vcf.gz', '')

filename_out += '+NoMasked.vcf'
f_out = open(filename_out, 'w')

count_var = 0

count_masked = 0
count_not_masked = 0
for line in f_vcf:
    if line.startswith('#'):
        f_out.write(line.strip() + "\n")
        continue

    tokens = line.strip().split("\t")
    tmp_chr = tokens[0]
    tmp_chr_pos = int(tokens[1])

    if tmp_chr not in seq_dict:
        sys.stderr.write('Not available in masked fa: %s\n' % tmp_chr)
    else:
        if seq_dict[tmp_chr][tmp_chr_pos] != 'N':
            f_out.write(line.strip()+"\n")
            count_not_masked += 1
        else:
            count_masked += 1
f_vcf.close()
sys.stderr.write('Total: %d, Masked: %d, NotMasked: %d\n' %
                 (count_masked+count_not_masked,
                  count_masked, count_not_masked))
