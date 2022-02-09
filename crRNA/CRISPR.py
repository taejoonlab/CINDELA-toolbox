import sys
import re

class CRISPR:
    def __init__(self):
        self.length_with_PAM = 0
        self.length_without_PAM = 0
        self.max_mismatch = 1

    def revcomp(self, tmp_seq):
        rc = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return ''.join([rc[x] for x in tmp_seq[::-1]])

    def validate_crRNA(self, tmp_seq):
        rv = False

        tmp_crRNA = tmp_seq.replace(' ', '')
        tmp_crRNA = tmp_crRNA.upper()
        if len(set(list(tmp_crRNA)) - set(['A', 'T', 'G', 'C'])) > 0:
            return False

        tmp_GC = tmp_seq.count('G') + tmp_crRNA.count('C')
        tmp_GC_ratio = tmp_GC * 100.0 / self.length_with_PAM
        if tmp_GC_ratio > 60 and tmp_GC_ratio < 40:
            return False
        return True
    
    def has_PAM(self, tmp_seq):
        return False

    def has_rcPAM(self, tmp_seq):
        return False

    def check_PAM(self, tmp_seq, reverse=False):
        tmp_seq = tmp_seq.replace(' ', '')
        if reverse is False:
            return has_PAM(tmp_seq)
        else:
            return has_rcPAM(tmp_seq)
        return False
    
    def max_homopolymer(self, tmp_seq):
        max_homopolymer = 1
        for tmp_n in ['A', 'T', 'G', 'C']:
            regex_n = re.compile('[%s]+' % tmp_n)
            for tmp_nx in regex_n.findall(tmp_seq.upper()):
                if len(tmp_nx) > max_homopolymer:
                    max_homopolymer = len(tmp_nx)
        return max_homopolymer

    def read_exonerate(self, tmp_filename):
        rv = dict()
        sys.stderr.write('Read %s (mm=%d)\n' % (tmp_filename, self.max_mismatch))
        f_exo_out = open(tmp_filename, 'r')
        for line in f_exo_out:
            if line.startswith('vulgar:'):
                tokens = line.strip().split()
                tmp_q_id = tokens[1]
                tmp_match_len = int(tokens[11])
                if tmp_match_len < self.length_with_PAM - self.max_mismatch:
                    continue
                if tmp_q_id not in rv:
                    rv[tmp_q_id] = []
                # tokens[7] : target_id
                rv[tmp_q_id].append(tokens[7])
        f_exo_out.close()
        return rv


class SpCas9(CRISPR):
    def __init__(self):
        super().__init__()
        # the length with PAM site: 20 + 3
        self.length_with_PAM = 23
        self.length_without_PAM = 20
        self.version = 'SpCRv4'
    
    def has_PAM(self, tmp_seq):
        tmp_pam = tmp_seq[-3:].upper()
        if tmp_pam.endswith('GG'):
            return True
        return False

    def has_rcPAM(self, tmp_seq):
        tmp_rc_pam = tmp_seq[:3].upper()
        if tmp_rc_pam.startswith('CC'):
            return True
        return False


class SaCas9(CRISPR):
    def __init__(self):
        super().__init__()
        # the length with PAM site: 21 + 6
        self.length_with_PAM = 27
        self.length_without_PAM = 21
        self.version = 'SaCRv4'

    def has_PAM(self, tmp_seq):
        # NNGRRT
        tmp_pam = tmp_seq[-6:].upper()
        if tmp_pam[2] == 'G' and tmp_pam[5] == 'T':
            if len(set(tmp_pam[3:5]) - set(['A','G'])) == 0:
                return True
        return False
   
    def has_rcPAM(self, tmp_seq):
        # AYYCNN
        tmp_rc_pam = tmp_seq[:6].upper()
        if tmp_rc_pam[0] == 'A' and tmp_rc_pam[3] == 'C':
            if len(set(tmp_rc_pam[1:3]) - set(['T','C'])) == 0:
                return True
        return False
