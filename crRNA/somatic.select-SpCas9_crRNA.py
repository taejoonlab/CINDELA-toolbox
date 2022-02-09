#!/usr/bin/env python3
import sys

len_crRNA = 23
max_mismatch = 2

# from var_flank_seq-to-crRNA.py
filename_cr_fa = sys.argv[1]

# by running exonerate
filename_exo_ref = filename_cr_fa.replace('.fa', '.REF.vulgar')
filename_exo_germline = filename_cr_fa.replace('.fa', '.GERMLINE.vulgar')
filename_exo_somatic = filename_cr_fa.replace('.fa', '.SOMATIC.vulgar')

chr_var_list = dict()
var_freq = dict()


def check_SpCas9_PAM(tmp_seq, reverse=False):
    tmp_seq = tmp_seq.replace(' ','')
    is_pam = False 
    if reverse == False:
        tmp_pam = tmp_seq[-3:]
        # NGG
        if tmp_pam[1] == 'G' and tmp_pam[2] == 'G':
            is_pam = True
    else: 
        # CCN
        tmp_rc_pam = tmp_seq[:3]
        if tmp_rc_pam[0] == 'C' and tmp_rc_pam[1] == 'C':
            is_pam = True

    return is_pam 


var2crRNA = dict()
var2score = dict()
crRNA_list = dict()

#>SpCRv2.AMC_CRC_BMK.001|chr1:842985:AACAGCTGAG-A|chr1_842K:AACAGCTGAG-A|F|18|TIR=45_0:10_34
#AGGATGGGAGAAGGTGATACTGG

is_select = 0
f_fa = open(filename_cr_fa, 'r')
for line in f_fa:
    if line.startswith('>'):
        seq_h = line.strip().lstrip('>')

        tmp_var_id = seq_h.split('|')[1]
        tmp_var_pos = ':'.join(tmp_var_id.split(':')[:2])
        tmp_var_TIR = seq_h.split('|')[-1].split('=')[1]
        tmp_ref_tokens = tmp_var_TIR.split(':')[0].split('_')
        tmp_alt_tokens = tmp_var_TIR.split(':')[1].split('_')
        tmp_var_score = int(tmp_alt_tokens[1]) - int(tmp_alt_tokens[0])

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


def read_exonerate(tmp_filename, tmp_mismatch):
    rv = []
    sys.stderr.write('Read %s (mm=%d)\n' % (tmp_filename, tmp_mismatch))
    f_exo_out = open(tmp_filename, 'r')
    for line in f_exo_out:
        if line.startswith('vulgar:'):
            tokens = line.strip().split()
            tmp_q_id = tokens[1]
            tmp_match_len = int(tokens[11])
            if tmp_match_len < len_crRNA - tmp_mismatch:
                continue
            rv.append(tmp_q_id)
    f_exo_out.close()
    return list(set(rv))


matched_var_ref_list      = read_exonerate(filename_exo_ref, max_mismatch)
matched_var_germline_list = read_exonerate(filename_exo_germline, 0)
matched_var_somatic_list  = read_exonerate(filename_exo_somatic, 0)


rc = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
def revcomp(tmp_seq):
    return ''.join([rc[x] for x in tmp_seq[::-1]])


for tmp_var_id in sorted(var2score.keys(), key=var2score.get, reverse=True):
#for tmp_var_id in sorted(var2crRNA.keys()):
    max_pos = -1
    max_pos_h = ''
    for tmp_h in var2crRNA[tmp_var_id]:
        if tmp_h in matched_var_ref_list:
            sys.stderr.write('%s in REF.fa\n' % tmp_h)
            continue
        
        if tmp_h in matched_var_germline_list:
            sys.stderr.write('%s in GERMLINE.fa\n' % tmp_h)
            continue

        count_somatic_hit = matched_var_somatic_list.count(tmp_h)
        if count_somatic_hit < 1:
            continue

        tmp_pos = int(tmp_h.split('|')[-2])
        if tmp_pos > max_pos:
            max_pos = tmp_pos
            max_pos_h = tmp_h
    
    if max_pos < 0:
        #sys.stderr.write('Failed: %s\n' % tmp_var_id)
        continue

    tmp_seq = crRNA_list[max_pos_h]

    if check_SpCas9_PAM(tmp_seq, reverse=False):
        print(">%s\n%s" % (max_pos_h, tmp_seq))
        #print(">%s\n%s" % (max_pos_h.split('|')[0], tmp_seq[:20]))

    elif check_SpCas9_PAM(tmp_seq, reverse=True):
        print(">%s revcomp\n%s" % (max_pos_h, revcomp(tmp_seq)))
        #print(">%s\n%s" % (max_pos_h.split('|')[0], revcomp(tmp_seq)[:20]))

    else:
        sys.stderr.write('PAM error: %s %s\n' % (max_pos_h, tmp_seq))
