#!/usr/bin/env python3
import sys
from CRISPR import SpCas9, SaCas9

usage = '%s <var_flank_seq> <target_name> <SpCas9 or SaCas9>' % sys.argv[0]

if len(sys.argv) != 4:
    sys.stderr.write('Usage: %s\n' % usage)
    sys.exit()

filename_var_flank_seq = sys.argv[1]
target_name = sys.argv[2]
crRNA_type = sys.argv[3]

if crRNA_type == 'SpCas9':
    cr = SpCas9()
elif crRNA_type == 'SaCas9':
    cr = SaCas9()

len_crRNA = cr.length_with_PAM

crRNA_idx = 1
f_var_flank = open(filename_var_flank_seq, 'r')
for line in f_var_flank:
    if line.startswith('#'):
        continue
    tokens = line.strip().split("\t")
    tmp_chr = tokens[0]
    tmp_pos = int(tokens[1])
    tmp_ref_allele = tokens[2]
    tmp_alt_allele = tokens[3]

    tmp_AD = tokens[4]
    tmp_ref_seq = tokens[5].replace('.', '')
    tmp_alt_seq_orig = tokens[6]
    var_pos_orig = tmp_alt_seq_orig.find('.')
    tmp_alt_seq = tmp_alt_seq_orig.replace('.', '')

    tmp_var_id = '%s:%d:%s-%s' %\
                 (tmp_chr, tmp_pos, tmp_ref_allele, tmp_alt_allele)

    for i in range(0, len(tmp_alt_seq)-len_crRNA+1):
        tmp_crRNA = tmp_alt_seq[i:i+len_crRNA]

        if cr.validate_crRNA(tmp_crRNA) is False:
            continue

        tmp_var_start = -1

        if cr.has_PAM(tmp_crRNA):
            tmp_var_start = var_pos_orig - i + 1
            #tmp_var_end = tmp_var_start + len(tmp_alt_allele)
            tmp_h_crRNA = '%s.%s.%03d-F|%s' %\
                          (cr.version, target_name, crRNA_idx, tmp_var_id)
            crRNA_idx += 1

        if cr.has_rcPAM(tmp_crRNA):
            tmp_var_start = i + len_crRNA - var_pos_orig - len(tmp_alt_allele) - 1
            #tmp_var_end = i + len_crRNA - var_pos_orig
            tmp_h_crRNA = '%s.%s.%03d-R|%s' %\
                          (cr.version, target_name, crRNA_idx, tmp_var_id)
            crRNA_idx += 1

        if tmp_ref_seq[i:i+len_crRNA] == tmp_crRNA:
            continue
        if tmp_var_start < 6 or tmp_var_start > len_crRNA:
            continue
        
        print(tmp_alt_seq_orig)
        print(tmp_var_start, var_pos_orig, i, tmp_crRNA, tmp_ref_seq[i:i+len_crRNA])

        if tmp_crRNA.find('AAAAAAA') >= 0 or tmp_crRNA.find('TTTTTTT') >= 0 \
            or tmp_crRNA.find('GGGGGGG') >= 0 or tmp_crRNA.find('CCCCCCC') >= 0:
            continue

        print(">%s|%d|AD=%s\n%s" %
              (tmp_h_crRNA, tmp_var_start, tmp_AD, tmp_crRNA))

f_var_flank.close()
