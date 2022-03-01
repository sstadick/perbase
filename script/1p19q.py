# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 14:04:42 2021

@author: xijun
"""


import sys
import numpy as np


def bed2pos(rs_pos_bed):
    pos_dict = {}
    pos_list = []
    with open(rs_pos_bed) as f_in:
        data_list = [i.split('\t') for i in f_in.read().strip().split('\n')]
        for line_list in data_list:
            pos = line_list[2]
            pos_list.append(pos)
            pos_dict[pos] = line_list
    return pos_list, pos_dict


def generate_writeline(data_perbase_list, rs_pos_list, rs_pos_dict):
    out_list = [['CHR', 'POS', 'DEPTH', 'REF', 'ALT', 'REF_DEPTH', 'ALT_DEPTH', 'FREQ', 'GENOTYPE']]
    out_dict = {}
    for pos in rs_pos_list:
        REF, pos_start, POS, REF_BASE = rs_pos_dict[pos]
        new_line_list = [REF, POS, '0', REF_BASE, REF_BASE, '0', '0', 'NA', '2']
        out_dict[pos] = new_line_list
        for line_list in data_perbase_list:
            if line_list[0] not in  ['1', '19']:
                continue
            if pos in line_list:
                REF, POS, REF_BASE, DEPTH, A, C, G, T, N = line_list[:9]
                new_line_list = [REF, POS, DEPTH, REF_BASE, REF_BASE, '0', '0', 'NA', '2']

                REF, POS, DEPTH, A, C, G, T, N = [int(i) for i in [REF, POS, DEPTH, A, C, G, T, N]]
                base_depth_list = [A, C, G, T, N]
                base_name_list = ['A', 'C', 'G', 'T', 'N']
                ref_depth = base_depth_list[base_name_list.index(REF_BASE)]

                max_depth = max(base_depth_list)
                max_depth_base = base_name_list[np.argmax(base_depth_list)]

                max2_depth = sorted(base_depth_list)[-2]
                max2_depth_base = base_name_list[base_depth_list.index(max2_depth)]

                # find alt and alt_depth
                if REF_BASE == max_depth_base:
                    alt = max2_depth_base
                    alt_depth = max2_depth
                    if max2_depth == max_depth:
                        alt = 'N'
                        if max_depth_base == REF:
                            alt = max2_depth_base
                        elif max2_depth_base == REF:
                            alt = max_depth_base
                else:  # germline
                    alt = max_depth_base
                    alt_depth = max_depth

                # calculate freq
                freq = (alt_depth/DEPTH)
                if DEPTH == 0:
                    freq = '0'
                # judge genotype
                if freq < 0.05 or freq > 0.95:
                    genotype = '2'
                elif 0.40 < freq < 0.60:
                    genotype = '0'
                else:
                    genotype = '1'
                new_line_list = [REF, POS, DEPTH, REF_BASE, alt, ref_depth, alt_depth, freq, genotype]
        out_list.append(new_line_list)
    return out_list


def write_list_out(sta_list, out_file):
    with open(out_file, 'w') as f_out:
        for line_list in sta_list:
            out_line = '\t'.join([str(i) for i in line_list])
            f_out.write(out_line + '\n')


if __name__ == '__main__':
    rs_pos_bed = sys.argv[1]  # '1p19q.rs.bed'
    perbase_file = sys.argv[2]  # 'out'
    out_file = sys.argv[3]  # 'N_1p19q.result.xls'

    rs_pos_list, rs_pos_dict = bed2pos(rs_pos_bed)

    out_list = []
    with open(perbase_file, 'r') as f_in:
        data_perbase_list = [i.split('\t') for i in f_in.read().strip().split('\n')]
    out_list = generate_writeline(data_perbase_list, rs_pos_list, rs_pos_dict)
    write_list_out(out_list, out_file)

