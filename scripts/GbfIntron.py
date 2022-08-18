#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
import os
import sys
import re

# Extract intron from GenBank file
# intron_gff.py | sort -k 1,1 -k 4n,4 -k 5n,5 | uniq | cut -f 1,4,5 > outfile
infile = sys.argv[1]

with open(infile, mode='r') as fi:
    entries = fi.read().split('//\n')
    for entry in entries:
        if len(entry) > 100:
            l_entry = entry.split('/')
            acc = re.search('LOCUS *.*?\n', entry).group().rstrip().split(' ')[7]   ####    acc without ver.
            # acc = re.search('\nVERSION *.*?\n', entry).group().rstrip().split(' ')[-1]    ####    acc with ver.
            for i in l_entry:
                if '   CDS   ' in i:
                    if 'complement(' in i:
                        strand = '-'
                    else:
                        strand = '+'
                    exons = re.sub(' |\n|complement|join|\(|\)', '', i.split('   CDS   ')[1])
                    l_exons = exons.split(',')

                    #### intron block
                    if len(l_exons) > 1:
                        for j in range(len(l_exons)):
                            if j == 0:
                                intron_start = int(re.sub('<|>', '', l_exons[j].split('..')[-1])) + 1
                            else:
                                intron_end = int(re.sub('<|>', '', l_exons[j].split('..')[0])) - 1
                                print('{}\tGenBank\tintron\t{}\t{}\t.\t{}\t.'.format(acc, intron_start, intron_end, strand))# output intron
                                intron_start = int(re.sub('<|>', '', l_exons[j].split('..')[-1])) + 1

