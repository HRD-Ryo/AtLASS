#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import subprocess
import sys
import random


intron = sys.argv[1]
exon = sys.argv[2]
posi = sys.argv[3]
nega = sys.argv[4]
min_intron = 20
min_exon = 20


with open(intron, mode='r') as fi, open(posi, mode='w') as fp, open(nega, mode='w') as fn:
    for i in fi:
        l_line = i.rstrip().split('\t')
        intron_len = int(l_line[2]) - int(l_line[1])
        if intron_len > min_intron:
            fp.write('{}\t{}\tforward\n'.format(l_line[0], l_line[1]))
            fn.write('{}\t{}\trevcom\n'.format(l_line[0], l_line[1]))
            fp.write('{}\t{}\trevcom\n'.format(l_line[0], l_line[2]))
            fn.write('{}\t{}\tforward\n'.format(l_line[0], l_line[2]))
            c = int(intron_len/100) + 1    ####    cycle num
            for j in range(c):
                p = random.randrange(int(l_line[1])+1, int(l_line[2]))
                fn.write('{}\t{}\tforward\n'.format(l_line[0], p))
                p = random.randrange(int(l_line[1])+1, int(l_line[2]))
                fn.write('{}\t{}\trevcom\n'.format(l_line[0], p))

with open(exon, mode='r') as fe, open(nega, mode='a') as fn:
    for i in fe:
        l_line = i.rstrip().split('\t')
        exon_len = int(l_line[2]) - int(l_line[1])
        if exon_len > min_exon:
            c = int(exon_len/100) + 1    ####    cycle num
            for j in range(c):
                p = random.randrange(int(l_line[1]), int(l_line[2])+1)
                fn.write('{}\t{}\tforward\n'.format(l_line[0], p))
                p = random.randrange(int(l_line[1]), int(l_line[2])+1)
                fn.write('{}\t{}\trevcom\n'.format(l_line[0], p))

