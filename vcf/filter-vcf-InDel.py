#!/usr/bin/env python3
import gzip
import sys

min_indel_len = 1
max_indel_len = 9

filename_vcf = sys.argv[1]
filename_out = filename_vcf.replace('.vcf', '')

f_vcf = open(filename_vcf, 'r')
if filename_vcf.endswith('.gz'):
    f_vcf = gzip.open(filename_vcf, 'rt')
    filename_out = filename_vcf.replace('.vcf.gz', '')

count_total = 0
count_pass = 0
count_selected = 0

filename_out += ('.+InDel%d-%d.vcf' % (min_indel_len, max_indel_len))
f_out = open(filename_out, 'w')

for line in f_vcf:
    if line.startswith('#'):
        f_out.write(line.strip() + "\n")
        continue

    count_total += 1

    tokens = line.strip().split("\t")
    tmp_ref = tokens[3]
    tmp_alt = tokens[4]

    # Do not count the long variants
    #if len(tmp_ref) > max_indel_len or max(tmp_len_alt) > max_indel_len:
    #    continue

    # Do not count on InDels with no call.
    if tmp_ref.find('N') >= 0 or tmp_alt.find('N') >= 0:
        continue

    # Do not count on InDels with ambiguity.
    if tmp_ref.find(',') >= 0 or tmp_alt.find(',') >= 0:
        continue

    #if tmp_pass_type.find('PASS') >= 0:
    #    count_pass += 1
    #else:
    #    continue
    
    tmp_pass_type = tokens[6]

    tmp_len_diff = abs(len(tmp_ref) - len(tmp_alt))
    if tmp_len_diff > max_indel_len or len(tmp_alt) > max_indel_len:
        continue

    # Do not count a InDel which is too long/short.
    if tmp_len_diff >= min_indel_len and tmp_len_diff <= max_indel_len:
        #tmp_ref_n = len(set(tmp_ref))
        #tmp_alt_n = max([len(set(tmp_seq)) for tmp_seq in tmp_alt.split(',')])
        # Skip homopolymer candidates
        # Removed. It cannot capture the homopolymer. - 2021.05.24.
        #if max([tmp_ref_n, tmp_alt_n]) < 4:
        #    continue

        count_selected += 1
        f_out.write(line.strip() + "\n")
f_out.close()

sys.stderr.write('Total variants: %d\n' % count_total)
sys.stderr.write('Passed variants: %d\n' % count_pass)
sys.stderr.write('Filtered variants: %d\n' % count_selected)
