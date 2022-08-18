#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
import os
import sys
import fileinput


def cigar(cigar, min_match, max_intron):
    length = ''
    M_length = 'ok'
    N_length = 'ok'
    for i in list(cigar):
        if i.isdecimal():
            length = length + i
        elif i == 'M':
            if int(length) < min_match:
                M_length = 'ng'
            length = ''
        elif i == 'N':
            if int(length) > max_intron:
                N_length = 'ng'
            length = ''
        else:
            length = ''
    return M_length, N_length


def ExtractReads():
    flag_normal = [99, 147, 83, 163]
    min_match = 10  ####    min length of match
    max_intron = 5000   ####    max size of cis-intron
    min_mapq = 3    ####    min MAPq
    pre_line = ""
    pre_l_list = [""]
    ####    print header
    for line in sys.stdin:
        if line.startswith('@'):
            print(line.rstrip())
    ####    print reads
        elif line != '':
            l_list = line.rstrip().split('\t')
            M_length, N_length = cigar(l_list[5], min_match, max_intron)
            if int(l_list[4]) >= min_mapq and M_length == 'ok' and N_length == 'ok' and int(l_list[1]) in flag_normal: #### MAPQc & cigar filter & flag
                print(line.rstrip())


ExtractReads()