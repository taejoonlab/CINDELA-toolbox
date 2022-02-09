#!/usr/bin/env python3
import sys
import gzip

filename_flank_seq = sys.argv[1]

f_flank_seq = open(filename_flank_seq, 'r')
if filename_flank_seq.endswith('.gz'):
    f_flank_seq = gzip.open(filename_flank_seq, 'rt')

filename_tokens = filename_flank_seq.split('.')
filename_out = filename_flank_seq.replace('.var_flank_seq', '') + ".alt.fa"

f_out = open(filename_out, 'w')
for line in f_flank_seq:
    if line.startswith('#'):
        continue
    tokens = line.strip().split("\t")
    tmp_chr_id = tokens[0]
    tmp_pos = tokens[1]
    tmp_ref = tokens[2]
    tmp_alt = tokens[3]
    tmp_var_id = '%s:%s:%s-%s' % (tmp_chr_id, tmp_pos, tmp_ref, tmp_alt)

    tmp_alt_tokens = tokens[6].split('.')
    tmp_alt_seq = ''.join(tmp_alt_tokens)
    if tmp_alt_seq.find('N') >= 0:
        continue

    f_out.write(">%s\n%s\n" % (tmp_var_id, tmp_alt_seq))
f_flank_seq.close()
f_out.close()
