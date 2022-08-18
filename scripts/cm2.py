#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
import os
import sys
import numpy as np
import math
from Bio import SeqIO
from sklearn.metrics import confusion_matrix


#summary
#   TN  FN  FP  TN  F1  MCC
# GMES  1780    876 502 1479023 0.720939651680842   0.722559030036
# GMET  1785    871 473 1479052 0.726495726495726   0.728445675655856
# fun_ES    1596    1060    360 1479165 0.692107545533391   0.699770806534218
# fun_ET    1759    897 399 1479126 0.730785209804736   0.734304387051734


genome_file = sys.argv[1]
true_file = sys.argv[2]
my_tsv_1 = sys.argv[3]
my_tsv_2 = sys.argv[4]

th1 = 0.5
th2 = 0.9

def f1_MCC(tp, fn, fp, tn):
    if tp == 0:
        f1, MCC = 0,0
    else:
        f1 = tp*2 / (2*tp+fp+fn)
        MCC = (tp*tn - fp*fn) / math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    return f1, MCC


d_genome_len = {}
for genome in SeqIO.parse(genome_file, 'fasta'):
    d_genome_len[genome.id] = len(genome.seq)


print('th1={}'.format(th1))
print('th2={}'.format(th2))
tp_all, fn_all, fp_all, tn_all = 0,0,0,0
for i in d_genome_len:
    l_true = np.zeros(d_genome_len[i])
    l_pred_1 = np.zeros(d_genome_len[i])
    l_pred_2 = np.zeros(d_genome_len[i])
    with open(true_file, mode='r') as f_true:
        for line in f_true:
            l_line = line.rstrip().split('\t')
            if l_line[0] == i:
                l_true[int(l_line[1])-1] = 1
                l_true[int(l_line[2])-1] = 1
    with open(my_tsv_1, mode='r') as f_my1:
        for line2 in f_my1:
            l_line2 = line2.rstrip().split('\t')
            if l_line2[0] == i and (float(l_line2[2]) > th1 or float(l_line2[3]) > th1):
                l_pred_1[int(l_line2[1])-1] = 1
    with open(my_tsv_2, mode='r') as f_my2:
        for line2 in f_my2:
            l_line2 = line2.rstrip().split('\t')
            if l_line2[0] == i and (float(l_line2[2]) > th2 or float(l_line2[3]) > th2):
                l_pred_2[int(l_line2[1])-1] = 1

    l_pred = (l_pred_1 * l_pred_2)


    cm = confusion_matrix(l_true, l_pred, labels=[1,0])
    # print(cm.ravel())
    tp, fn, fp, tn = cm.ravel()
    tp_all += tp
    fn_all += fn
    fp_all += fp
    tn_all += tn
print(tp_all, fn_all, fp_all, tn_all)
print(f1_MCC(int(tp_all), int(fn_all), int(fp_all), int(tn_all)))


# .as(int)