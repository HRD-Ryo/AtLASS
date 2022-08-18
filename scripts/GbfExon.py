#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
import os
import sys
import re

# Extract Exon from GenBank file
# intron_gff.py | sort -k 1,1 -k 4n,4 | uniq > outfile
infile = sys.argv[1]

with open(infile, mode='r') as fi:
    entries = fi.read().split('//\n')
    for entry in entries:
        if entry != '':
            l_entry = entry.split('/')
            # acc = re.search('LOCUS *.*?\n', entry).group().rstrip().split(' ')[7]   ####    acc without ver.
            acc = re.search('\nVERSION *.*?\n', entry).group().rstrip().split(' ')[-1]    ####    acc with ver.
            for i in l_entry:
                if '   CDS   ' in i:
                    if 'complement(' in i:
                        strand = '-'
                    else:
                        strand = '+'
                    exons = re.sub(' |\n|complement|join|\(|\)', '', i.split('   CDS   ')[1])
                    l_exons = exons.split(',')

                    #### exon block
                    for k in l_exons:
                        l_exons_2 = k.split('..')
                        if len(l_exons_2) == 2:
                            print('{}\t{}\t{}'.format(acc, re.sub('<|>', '', l_exons_2[0]), re.sub('<|>', '', l_exons_2[1])))

                    #### intron block
                    # if len(l_exons) > 1:
                    #     for j in range(len(l_exons)):
                    #         if j == 0:
                    #             intron_start = int(re.sub('<|>', '', l_exons[j].split('..')[-1])) + 1
                    #         else:
                    #             intron_end = int(re.sub('<|>', '', l_exons[j].split('..')[0])) - 1
                    #             print('{}\tGenBank\tintron\t{}\t{}\t.\t{}\t.'.format(acc, intron_start, intron_end, strand))# output intron
                    #             intron_start = int(re.sub('<|>', '', l_exons[j].split('..')[-1])) + 1


