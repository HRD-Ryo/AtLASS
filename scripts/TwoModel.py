#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import subprocess
import sys
import pandas as pd
from Bio import SeqIO


model1 = sys.argv[1]
model2 = sys.argv[2]
out = sys.argv[3]
tmp1 = '/tmp/TwoModel_tmp1'
tmp2 = '/tmp/TwoModel_tmp2'

th1 = 0.5
th2 = 0.5

def Include(f, df, sca, posi, score1, score2):
    if score1 > th1:
        df_tmp = df[(df[0]==sca) & (df[1]==posi) & (df[2]>th2)]
        if len(df_tmp) > 0:
            f.write('{}\t{}\t{}\t{}\n'.format(sca, posi, 'forward', round((score1+df_tmp.iat[0,2])/2, 3)))
    if score2 > th1:
        df_tmp = df[(df[0]==sca) & (df[1]==posi) & (df[3]>th2)]
        if len(df_tmp) > 0:
            f.write('{}\t{}\t{}\t{}\n'.format(sca, posi, 'revcom', round((score2+df_tmp.iat[0,3])/2, 3)))


df = pd.read_csv(model2, header=None, sep='\t')
with open(tmp1, mode='w') as ft1, open(model1, mode='r') as f1:
    for i in f1:
        line = i.rstrip().split('\t')
        sca, posi, score1, score2 = line[0], int(line[1]), float(line[2]), float(line[3])
        Include(ft1, df, sca, posi, score1, score2)


pre = ['', 0, '', 0]
with open(tmp2, mode='w') as ft2, open(tmp1, mode='r') as ft1:
    for i in ft1:
        line = i.rstrip().split('\t')
        if line[0] == pre[0] and line[2] == pre[2] and int(line[1]) - int(pre[1]) < 2:
            if float(line[3]) > float(pre[3]):
                ft2.write('\t'.join(pre)+'\n')
            else:
                ft2.write('\t'.join(line)+'\n')
        pre = line


c = 'grep -v -f {} {} > {}'.format(tmp2, tmp1, out)
os.system(c)
os.remove(tmp1)
os.remove(tmp2)
