#!/usr/bin/env python3
import sys
import gzip

filename_vcf = sys.argv[1]

# Make a symbolic link for the reference FASTA to 'REF.fa'
filename_ref_fa = "REF.fa"

# Flanking sequences from the InDel position (both direction)
len_flank = 30

refseq_list = dict()
f_ref_fa = open(filename_ref_fa, 'r')
if filename_ref_fa.endswith('.gz'):
    f_ref_fa = gzip.open(filename_ref_fa, 'rt')

for line in f_ref_fa:
    if line.startswith('>'):
        seq_h = line.strip().split()[0].lstrip('>')
        refseq_list[seq_h] = []
    else:
        refseq_list[seq_h].append(line.strip())
f_ref_fa.close()

refseq_concat = dict()
for tmp_h in refseq_list.keys():
    refseq_concat[tmp_h] = ''.join(refseq_list[tmp_h])

f_vcf = open(filename_vcf, 'r')
if filename_vcf.endswith('.gz'):
    f_vcf = gzip.open(filename_vcf, 'rt')

print("#CHROM\tPOS\tREF\tALT\tDP\tREF_flank\tALT_flank")
for line in f_vcf:
    if line.startswith('#'):
        continue
    tokens = line.strip().split("\t")
    tmp_seq_id = tokens[0]
    tmp_seq_pos = int(tokens[1])
    tmp_ref_allele = tokens[3]
    tmp_alt_allele = tokens[4]
    if tmp_alt_allele.find(',') >= 0:
        continue

    fmt_tokens = tokens[8].split(':')
    idx_DP = fmt_tokens.index('DP')

    normal_fmt_tokens = tokens[9].split(':')
    tumor_fmt_tokens = tokens[10].split(':')
    normal_DP = int(normal_fmt_tokens[idx_DP])
    tumor_DP = int(tumor_fmt_tokens[idx_DP])

    tmp_DP= '%d:%d' % (normal_DP, tumor_DP)

    len_ref_allele = len(tmp_ref_allele)
    tmp_len_flank = len_flank
    tmp_ref = refseq_concat[tmp_seq_id][tmp_seq_pos-len_flank:tmp_seq_pos+len_flank]
    tmp_ref_5 = refseq_concat[tmp_seq_id][tmp_seq_pos-len_flank:tmp_seq_pos-1]
    tmp_ref_0 = refseq_concat[tmp_seq_id][tmp_seq_pos-1:tmp_seq_pos+len_ref_allele-1]
    tmp_ref_3 = refseq_concat[tmp_seq_id][tmp_seq_pos+len_ref_allele-1:tmp_seq_pos+len_flank]

    tmp_ref_seq = '%s.%s.%s' % (tmp_ref_5, tmp_ref_allele, tmp_ref_3)
    tmp_alt_seq = '%s.%s.%s' % (tmp_ref_5, tmp_alt_allele, tmp_ref_3)

    if tmp_ref_seq.count('N') > 10:
        continue
    
    print("%s\t%d\t%s\t%s\t%s\t%s\t%s" % (tmp_seq_id, tmp_seq_pos, tmp_ref_allele, tmp_alt_allele, tmp_DP, tmp_ref_seq, tmp_alt_seq))
f_vcf.close()
