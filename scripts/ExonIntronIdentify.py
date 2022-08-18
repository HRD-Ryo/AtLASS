#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
import os
import sys
import fileinput

####    input == stdin from mpileup
####    output1 == exon.tsv
####    output2 == exon.tsv
####    samtools mpileup -f genome bam | ExonIdentify.py exon.tsv intron.tsv

min_cov = 10
exon_tsv = sys.argv[1]
intron_tsv = sys.argv[2]

pre_scaffold = ''
pre_posi = 0
n_rate = 0
pre_n_rate = 0
stat_1 = 'normal'
stat_2 = 'normal'


with open(exon_tsv, mode='w') as fe, open(intron_tsv, mode='w') as fi:
    for line in sys.stdin:
            l_line = line.split('\t')
            if int(l_line[3]) > min_cov:
                if (l_line[4].count('>') + l_line[4].count('<')) == 0:
                    n_rate = 0
                else:
                    n_rate = (l_line[4].count('>') + l_line[4].count('<')) / int(l_line[3]) ####    rate of cigarN (gap)
                n_rate_dif = n_rate - pre_n_rate
                if pre_scaffold == l_line[0] and pre_posi == int(l_line[1])-1:  ####    continue

        ####    intron
                    if stat_2 != 'intron' and n_rate_dif > 0.5:   ####    start of intron
                        intron_start = l_line[1]
                        stat_2 = 'intron'
                        score1 = n_rate_dif
                    elif stat_2 == 'intron' and (n_rate_dif < -0.5 or n_rate < 0.5):  ####    end of intron
                        if score1 > 0.8 and -n_rate_dif > 0.8:
                            fi.write('{}\t{}\t{}\n'.format(pre_scaffold, intron_start, pre_posi))
                        stat_2 = 'normal'

        ####    exon
                    if stat_1 != 'exon' and n_rate < 0.2:
                        exon_start = l_line[1]
                        stat_1 = 'exon'
                    elif stat_1 == 'exon' and n_rate > 0.2:
                        stat_1 = 'normal'
                        fe.write('{}\t{}\t{}\n'.format(pre_scaffold, exon_start, pre_posi))
                elif stat_1 == 'exon':
                    stat_1 = 'normal'
                    fe.write('{}\t{}\t{}\n'.format(pre_scaffold, exon_start, pre_posi))

                else:
                    stat_1 = 'normal'
                    stat_2 = 'normal'

                pre_scaffold = l_line[0]
                pre_posi = int(l_line[1])
                pre_n_rate = n_rate