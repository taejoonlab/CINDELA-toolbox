#!/usr/bin/env python3
import sys

filename_var_flank_seq = sys.argv[1]
target_name = sys.argv[2]

len_crRNA= 23

#CHROM	POS	REF	ALT	TIR	REF_flank	ALT_flank
#chr1	597180	CTCAT	C	28_0:0_4	ATCCAAACACACAAAAATATATGCATGCG.CTCAT.TCATACACAATCTCACACATACATAT	ATCCAAACACACAAAAATATATGCATGCG.C.TCATACACAATCTCACACATACATAT

crRNA_idx = 1
f_var_flank = open(filename_var_flank_seq, 'r')
for line in f_var_flank:
    if line.startswith('#'):
        continue
    tokens = line.strip().split("\t")
    tmp_chr_id = tokens[0]
    if tmp_chr_id.startswith('chrUn'):
        continue

    tmp_pos = int(tokens[1])
    tmp_ref_allele = tokens[2]
    tmp_alt_allele = tokens[3]
    tmp_TIR = tokens[4]
    tmp_ref_seq = tokens[5].replace('.','')
    tmp_alt_seq_orig = tokens[6]
    tmp_alt_seq = tmp_alt_seq_orig.replace('.','')

    tmp_pos_tag = '%d' % tmp_pos
    if tmp_pos > 1000000:
        tmp_pos_tag = '%dM' % (tmp_pos / 1000000)
    elif tmp_pos > 1000:
        tmp_pos_tag = '%dK' % (tmp_pos / 1000)
    var_tag = '%s_%s:%s-%s' % (tmp_chr_id, tmp_pos_tag, tmp_ref_allele, tmp_alt_allele)

    var_id = 'SpCRv2.%s.%03d|%s:%d:%s-%s' % (target_name, crRNA_idx, tmp_chr_id, tmp_pos, tmp_ref_allele, tmp_alt_allele)
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
        
        ## PAM for SpCas9
        if tmp_crRNA.endswith('GG'):
            tmp_var_start = var_pos_orig - i
            tmp_var_end = tmp_var_start + len(tmp_alt_allele)
            if tmp_var_start < 6:
                continue
            print(">%s|%s|F|%d|TIR=%s\n%s" % (var_id, var_tag, tmp_var_start, tmp_TIR, tmp_crRNA))

        if tmp_crRNA.startswith('CC'):
            tmp_var_start = i + len_crRNA - var_pos_orig - len(tmp_alt_allele)
            tmp_var_end = i + len_crRNA - var_pos_orig
            if tmp_var_start < 6:
                continue
            print(">%s|%s|R|%d|TIR=%s\n%s" % (var_id, var_tag, tmp_var_start, tmp_TIR, tmp_crRNA))
f_var_flank.close()
