#!/usr/bin/env python3
import sys
import gzip

filename_vcf = sys.argv[1]

# Make a symbolic link for the reference FASTA to 'REF.fa'
filename_ref_fa = "REF.fa"

# Flanking sequences from the InDel position (both direction)
len_flank = 40

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

refseqs = dict()
for tmp_h in refseq_list.keys():
    refseqs[tmp_h] = ''.join(refseq_list[tmp_h])

f_vcf = open(filename_vcf, 'r')
if filename_vcf.endswith('.gz'):
    f_vcf = gzip.open(filename_vcf, 'rt')

var_id_list = []

print("#CHROM\tPOS\tREF\tALT\tAD\tREF_flank\tALT_flank")
for line in f_vcf:
    if line.startswith('#'):
        continue
    tokens = line.strip().split("\t")
    chr_id = tokens[0]
    tmp_pos = int(tokens[1])
    tmp_ref_allele = tokens[3]
    tmp_alt_allele = tokens[4]
    if tmp_alt_allele.find(',') >= 0:
        continue

    tmp_var_id = '%s|%d|%s-%s' % (chr_id, tmp_pos, tmp_ref_allele, tmp_alt_allele)
    if tmp_var_id in var_id_list:
        continue
    
    var_id_list.append(tmp_var_id)

    tmp_AD = 'NotAvail'
    #tmp_fmt_tokens = tokens[8].split(':')
    #idx_AD = tmp_fmt_tokens.index('AD')
    #tmp_AD = tokens[9].split(':')[idx_AD]

    len_ref_allele = len(tmp_ref_allele)
    tmp_ref = refseqs[chr_id][tmp_pos-len_flank:tmp_pos+len_flank]
    tmp_ref_5 = refseqs[chr_id][tmp_pos-len_flank-1:tmp_pos-1]
    tmp_ref_0 = refseqs[chr_id][tmp_pos-1:tmp_pos+len_ref_allele-1]
    tmp_ref_3 = refseqs[chr_id][tmp_pos+len_ref_allele-1:tmp_pos+len_ref_allele-1+len_flank]

    tmp_ref_seq = '%s.%s.%s' % (tmp_ref_5, tmp_ref_allele, tmp_ref_3)
    tmp_alt_seq = '%s.%s.%s' % (tmp_ref_5, tmp_alt_allele, tmp_ref_3)

    if tmp_ref_seq.count('N') > 10:
        continue

    print("%s\t%d\t%s\t%s\t%s\t%s\t%s" %
          (chr_id, tmp_pos, tmp_ref_allele, tmp_alt_allele,
           tmp_AD, tmp_ref_seq, tmp_alt_seq))
f_vcf.close()
