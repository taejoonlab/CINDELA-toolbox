#!/usr/bin/env python3
import sys
import gzip

filename_vcf = sys.argv[1]
min_AD = int(sys.argv[2])

prefix_out = '.baf.AD%d.bg' % min_AD
f_vcf = open(filename_vcf, 'r')
filename_out = filename_vcf.replace('.vcf', '') + prefix_out

if filename_vcf.endswith('.gz'):
    f_vcf = gzip.open(filename_vcf, 'rt')
    filename_out = filename_vcf.replace('.vcf.gz', '') + prefix_out

sys.stderr.write('Write %s\n' % filename_out)
f_out = open(filename_out, 'w')
for line in f_vcf:
    if line.startswith('#'):
        continue
    tokens = line.strip().split("\t")
    seq_id = tokens[0]
    tmp_pos = int(tokens[1])
    ref_allele = tokens[3]
    alt_allele = tokens[4]
    fmt_tokens = tokens[8].split(':')
    value_tokens = tokens[9].split(':')

    idx_AD = fmt_tokens.index('AD')
    tokens_AD = [int(x) for x in value_tokens[idx_AD].split(',')]

    sum_AD = sum(tokens_AD)
    if sum_AD < min_AD:
        continue

    tmp_BAF = tokens_AD[0] / sum_AD
    f_out.write("%s\t%d\t%d\t%.4f\n" %
                (seq_id, tmp_pos, tmp_pos+len(ref_allele), tmp_BAF))
f_vcf.close()
f_out.close()
