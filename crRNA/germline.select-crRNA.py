#!/usr/bin/env python3
import sys
from CRISPR import SpCas9, SaCas9

usage = '%s <FASTA for raw crRNAs> <SpCas9 or SaCas9>' % sys.argv[0]

if len(sys.argv) != 3:
    sys.stderr.write('Usage: %s\n' % usage)
    sys.exit()

# from var_flank_seq-to-crRNA.py
filename_cr_fa = sys.argv[1]
crRNA_type = sys.argv[2]

if crRNA_type == 'SpCas9':
    cr = SpCas9()
elif crRNA_type == 'SaCas9':
    cr = SaCas9()

len_crRNA = cr.length_with_PAM

# by running exonerate
filename_exo_ref = filename_cr_fa.replace('.fa', '.REF.vulgar')
filename_exo_germline = filename_cr_fa.replace('.fa', '.GERMLINE.vulgar')

chr_var_list = dict()
var_freq = dict()

var2crRNA = dict()
var2score = dict()
crRNA_list = dict()

is_select = 0
f_fa = open(filename_cr_fa, 'r')
for line in f_fa:
    if line.startswith('>'):
        seq_h = line.strip().lstrip('>')

        tmp_var_id = seq_h.split('|')[1]
        tmp_var_pos = ':'.join(tmp_var_id.split(':')[:2])
        tmp_var_AD = seq_h.split('|')[-1].split('=')[1]
        tmp_var_AD_ref = int(tmp_var_AD.split(',')[0])
        tmp_var_AD_germline = int(tmp_var_AD.split(',')[-1])
        tmp_var_score = tmp_var_AD_germline - tmp_var_AD_ref

        if tmp_var_pos.startswith('chrUn'):
            is_select = 0
        else:
            is_select = 1

        if is_select > 0:
            sys.stderr.write('Good: %s\n' % (tmp_var_pos))
            crRNA_list[seq_h] = ''
            if tmp_var_id not in var2crRNA:
                var2crRNA[tmp_var_id] = []
            var2crRNA[tmp_var_id].append(seq_h)
            var2score[tmp_var_id] = tmp_var_score

    elif is_select > 0:
        crRNA_list[seq_h] += line.strip()
f_fa.close()

matched_var_ref_list = cr.read_exonerate(filename_exo_ref)
matched_var_germline_list = cr.read_exonerate(filename_exo_germline)

count_total = 0
count_ref_match = 0
count_no_germline_match = 0
count_multi_germline_match = 0
count_selected = 0

for tmp_var_id in sorted(var2score.keys(), key=var2score.get, reverse=True):
    max_pos = -1
    max_pos_h = ''
    for tmp_h in var2crRNA[tmp_var_id]:
        count_total += 1
        if tmp_h in matched_var_ref_list:
            count_ref_match += 1
            sys.stderr.write('%s in REF.fa\n' % tmp_h)
            continue

        if tmp_h not in matched_var_germline_list:
            count_no_germline_match += 1
            continue
        else:
            count_germline_hit = len(set(matched_var_germline_list[tmp_h]))
            if count_germline_hit != 1:
                count_multi_germline_match += 1
                continue

        tmp_pos = int(tmp_h.split('|')[-2])
        if tmp_pos > max_pos and tmp_pos < len_crRNA - 3:
            max_pos = tmp_pos
            max_pos_h = tmp_h

    if max_pos < 0:
        sys.stderr.write('Failed: %s\n' % tmp_var_id)
        continue

    tmp_seq = crRNA_list[max_pos_h].upper()
    if cr.max_homopolymer(tmp_seq) > 5:
        sys.stderr.write('Max homopolymer: %s\n' % tmp_seq)
        continue

    if cr.has_PAM(tmp_seq):
        count_selected += 1
        print(">%s\n%s" % (max_pos_h, tmp_seq))

    elif cr.has_rcPAM(tmp_seq):
        count_selected += 1
        print(">%s revcomp\n%s" % (max_pos_h, cr.revcomp(tmp_seq)))
    else:
        sys.stderr.write('PAM error: %s %s\n' % (max_pos_h, tmp_seq))

sys.stderr.write('Total crRNAs: %d\n' % count_total)
sys.stderr.write('Matched to REF: %d\n' % count_ref_match)
sys.stderr.write('No matched to GERMLINE: %d\n' % count_no_germline_match)
sys.stderr.write('Multiple matched to GERMLINE: %d\n' %
                 count_multi_germline_match)
sys.stderr.write('Selected: %d\n' % count_selected)
