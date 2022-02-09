#!/usr/bin/env python3
import sys

filename_var_flank_seq = sys.argv[1]
target_name = sys.argv[2]

cr_version = 'SaCRv3'
len_crRNA= 27


def check_SaCas9_PAM(tmp_seq, reverse=False):
    tmp_seq = tmp_seq.replace(' ','')
    is_pam = False 
    if reverse == False:
        tmp_pam = tmp_seq[-6:]
        # NNGRRT
        if tmp_pam[2] == 'G' and tmp_pam[5] == 'T':
            if len(set(tmp_pam[3:5]) - set(['A','G'])) == 0:
                is_pam = True
    else: 
        # AYYCNN
        tmp_rc_pam = tmp_seq[:6]
        if tmp_rc_pam[0] == 'A' and tmp_rc_pam[3] == 'C':
            if len(set(tmp_rc_pam[1:3]) - set(['T','C'])) == 0:
                is_pam = True

    return is_pam 


#print(check_SaCas9_PAM('ATATGGATAT ATATAGAGGAGA GGGACT'))
#print(check_SaCas9_PAM('ATCCGG ATAT ATATAGAGGAGA GGGACT', reverse=True))
#sys.exit(1)

crRNA_idx = 1
f_var_flank = open(filename_var_flank_seq, 'r')
for line in f_var_flank:
    if line.startswith('#'):
        continue
    tokens = line.strip().split("\t")
    tmp_chr_id = tokens[0]
    tmp_pos = int(tokens[1])
    tmp_ref_allele = tokens[2]
    tmp_alt_allele = tokens[3]
    tmp_AD = tokens[4]
    tmp_ref_seq = tokens[5].replace('.','')
    tmp_alt_seq_orig = tokens[6]
    tmp_alt_seq = tmp_alt_seq_orig.replace('.','')
    var_id = '%s.%s.%03d|%s:%d:%s-%s' % (cr_version, target_name, crRNA_idx, tmp_chr_id, tmp_pos, tmp_ref_allele, tmp_alt_allele)
    crRNA_idx += 1
    
    var_pos_orig = tmp_alt_seq_orig.find('.')
    for i in range(0, len(tmp_alt_seq)-len_crRNA+1):
        tmp_crRNA = tmp_alt_seq[i:i+len_crRNA]
        if tmp_crRNA.find('N') >= 0:
            continue
        if tmp_ref_seq.find(tmp_crRNA) >= 0:
            continue
        
        tmp_GC = tmp_crRNA.count('G') + tmp_crRNA.count('C')
        tmp_GC_ratio = tmp_GC * 100.0 / len_crRNA
        if tmp_GC_ratio > 60 and tmp_GC_ratio < 40:
            continue
        
        ## PAM for SaCas9 : NNGRRT
        if check_SaCas9_PAM(tmp_crRNA, reverse=False):
            tmp_var_start = var_pos_orig - i
            tmp_var_end = tmp_var_start + len(tmp_alt_allele)
            if tmp_var_start < 6:
                continue
            print(">%s|F|%d|AD=%s\n%s" % (var_id, tmp_var_start, tmp_AD, tmp_crRNA))

        if check_SaCas9_PAM(tmp_crRNA, reverse=True):
            tmp_var_start = i + len_crRNA - var_pos_orig - len(tmp_alt_allele)
            tmp_var_end = i + len_crRNA - var_pos_orig
            if tmp_var_start < 6:
                continue
            print(">%s|R|%d|AD=%s\n%s" % (var_id, tmp_var_start, tmp_AD, tmp_crRNA))

f_var_flank.close()
