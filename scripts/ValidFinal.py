#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import subprocess
import sys
import pandas as pd
from Bio import SeqIO


####    validate GenomePred result
####    | sort -k 1,1 -k 2n,2


GenomePred = sys.argv[1]
genome_file = sys.argv[2]
intron = sys.argv[3]
exon = sys.argv[4]


def Include1(df_sca, posi, strand):
    df_tmp = df_sca[(df_sca[1]==posi) & (df_sca[2]==strand)]
    if len(df_tmp) > 0:
        r = 'True Positive'
    else:
        r = 'False Negative'
    return r


def Include2(df_intron, df_exon, posi, strand):
    if strand == 'forward':
        if len(df_intron[df_intron[1]==posi]) > 0:
            r = 'True Positive'
        elif len(df_intron[(df_intron[1]<posi) & (df_intron[2]>=posi)]) > 0 \
        or len(df_exon[(df_exon[1]<=posi) & (df_exon[2]>=posi)]) > 0:
            r = 'False Positive'
        else:
            r = 'AtLASS Prediction'
    elif strand == 'revcom':
        if len(df_intron[df_intron[2]==posi]) > 0:
            r = 'True Positive'
        elif len(df_intron[(df_intron[1]<=posi) & (df_intron[2]>posi)]) > 0 \
        or len(df_exon[(df_exon[1]<=posi) & (df_exon[2]>=posi)]) > 0:
            r = 'False Positive'
        else:
            r = 'AtLASS Prediction'
    return r


l_genome = []    ####    list of genome id
for genome in SeqIO.parse(genome_file, 'fasta'):
    l_genome.append(str(genome.id))

df_result = pd.read_csv(GenomePred, header=None, sep='\t')
for j in l_genome:
    df_sca = df_result[df_result[0]==j]
    with open(intron, mode='r') as f:
        for i in f:
            l_line = i.rstrip().split('\t')
            if l_line[0] == j:
                r = Include1(df_sca, int(l_line[1]), 'forward')
                if r == 'False Negative':
                    print('{}\t{}\tforward\t{}\t{}\n'.format(l_line[0], l_line[1], 'N/A', r), end='')
                r = Include1(df_sca, int(l_line[2]), 'revcom')
                if r == 'False Negative':
                    print('{}\t{}\trevcom\t{}\t{}\n'.format(l_line[0], l_line[2], 'N/A', r), end='')
del df_result


df_intron = pd.read_csv(intron, header=None, sep='\t')
df_exon = pd.read_csv(exon, header=None, sep='\t')
for j in l_genome:
    df_sca_intron = df_intron[df_intron[0]==j]
    df_sca_exon = df_exon[df_exon[0]==j]
    with open(GenomePred, mode='r') as f:
        for i in f:
            l_line = i.rstrip().split('\t')
            if l_line[0] == j:
                r = Include2(df_sca_intron, df_sca_exon, int(l_line[1]), l_line[2])
                print('{}\t{}\t{}\t{}\t{}\n'.format(l_line[0], l_line[1], l_line[2], l_line[3], r), end='')


