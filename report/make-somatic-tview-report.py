#!/usr/bin/env python3
import sys
import os

filename_fa = sys.argv[1]

rc = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}


def revcomp(tmp_seq):
    return ''.join([rc[x] for x in tmp_seq[::-1]])


def print_alignment(tmp_list):
    tmp_prev = ''
    tmp_count = 1
    for tmp_line in tmp_list:
        if tmp_line.count(' ') > len(tmp_line)*0.6:
            continue

        tmp_line = tmp_line.rstrip()
        if tmp_line == tmp_prev:
            tmp_count += 1
        else:
            print("%02d" % tmp_count, tmp_line)
            tmp_count = 1
        tmp_prev = tmp_line

    if tmp_line.count(' ') <= len(tmp_line)*0.6:
        print("%02d" % tmp_count, tmp_line)


f_fa = open(filename_fa, 'r')
for line in f_fa:
    if line.startswith('>'):
        tmp_cr_id = line.strip().lstrip('>')
        tmp_var_id = '%s.%s' % (tmp_cr_id.split('|')[0], tmp_cr_id.split('|')[1].replace(':', '_'))
        tmp_cr_seq = next(f_fa).strip()

        tview_list = {'NORMAL': [], 'TUMOR': []}
        for tmp_filename in os.listdir('.'):
            if tmp_filename.startswith(tmp_var_id):
                if tmp_filename.find('-normal') >= 0:
                    tmp_f = open(tmp_filename, 'r')
                    tview_list['NORMAL'] = tmp_f.readlines()
                    tmp_f.close()
                elif tmp_filename.find('-tumor') >= 0:
                    tmp_f = open(tmp_filename, 'r')
                    tview_list['TUMOR'] = tmp_f.readlines()
                    tmp_f.close()

        print("".join(["-=" for i in range(0, 40)]))
        print()

        print("crRNA ID: %s" % tmp_cr_id)
        print("crRNA seq: %s" % tmp_cr_seq)
        tmp_matched_seq = tmp_cr_seq
        if tmp_cr_id.find('revcomp') >= 0:
            tmp_matched_seq = revcomp(tmp_cr_seq)
        print("crRNA target matched seq: %s" % tmp_matched_seq)
        print()

        for tmp_type in ['NORMAL', 'TUMOR']:
            print("\n==== %s ====" % tmp_type)
            tmp_rep_seq = tview_list[tmp_type][1]
            max_pos = 0
            max_count_matched = 0
            # Position
            print(tview_list[tmp_type][0].strip())

            for i in range(0, len(tmp_rep_seq) - len(tmp_matched_seq)+1):
                count_matched = 0
                for j in range(0, len(tmp_matched_seq)):
                    if tmp_matched_seq[j] == tmp_rep_seq[i+j]:
                        count_matched += 1
                if count_matched > max_count_matched:
                    max_count_matched = count_matched
                    max_pos = i

            if count_matched > 15:
                print("   " + tmp_matched_seq.rjust(i+23, ' ') + " (matched: %d)" % count_matched)

            print("   " + tmp_matched_seq.rjust(max_pos+23, ' ') + " (best_matched: %d)" % max_count_matched)
            print_alignment(tview_list[tmp_type][1:])
f_fa.close()
