#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
import os
os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"    #### DDBJ
os.environ["CUDA_VISIBLE_DEVICES"]="0"    #### DDBJ
os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'max_split_size_mb:128'
import sys
from Bio import SeqIO
import numpy as np
import torch
import torch.nn as nn
import math
import gc
import re


genome_file = sys.argv[1]
model_path = sys.argv[2]
outfile = sys.argv[3]
if sys.argv[4] in ['True', 'TRUE', 'true']:
    gpu = True
else:
    gpu = False
device = torch.device('cuda' if torch.cuda.is_available() and gpu else 'cpu')
# print('{} is used'.format(device))

dim = 201
win = int((dim-1)/2)
####    fast mode
if sys.argv[5] in ['True', 'TRUE', 'true']:
    lstm_dim = 256
    lstm_layers = 1
    da = 128
    r = 1
    fc_dim = 32
    batch_size = 1000
else:
    lstm_dim = 1024
    lstm_layers = 1
    da = 1024
    r = 3
    fc_dim = 1024
    batch_size = 400


def onehot_plus(nucl, strand):
    if strand == 'forward':
        d_nucl = {'A':[1.,0.,0.,0.], 'T':[0.,1.,0.,0.], 'G':[0.,0.,1.,0.], 'C':[0.,0.,0.,1.]}
        all_encode = []
        for i in range(len(nucl)):
            s = nucl[i]
            onehot = d_nucl[s]
            freq_A = nucl[i:].count('A') / (len(nucl)-i)
            freq_T = nucl[i:].count('T') / (len(nucl)-i)
            freq_G = nucl[i:].count('G') / (len(nucl)-i)
            freq_C = nucl[i:].count('C') / (len(nucl)-i)
            onehot = onehot + [freq_A, freq_T, freq_G, freq_C]
            all_encode.append(onehot)
    elif strand == 'revcom':
        d_nucl = {'A':[0.,1.,0.,0.], 'T':[1.,0.,0.,0.], 'G':[0.,0.,0.,1.], 'C':[0.,0.,1.,0.]}
        nucl = nucl[::-1]
        all_encode = []
        for i in range(len(nucl)):
            s = nucl[i]
            onehot = d_nucl[s]
            freq_A = nucl[i:].count('T') / (len(nucl)-i)
            freq_T = nucl[i:].count('A') / (len(nucl)-i)
            freq_G = nucl[i:].count('C') / (len(nucl)-i)
            freq_C = nucl[i:].count('G') / (len(nucl)-i)
            onehot = onehot + [freq_A, freq_T, freq_G, freq_C]
            all_encode.append(onehot)
    return all_encode


class LSTM_Attn(nn.Module):
    def __init__(self, lstm_dim=lstm_dim, lstm_layers=lstm_layers, out_dim=2, da=da, r=r, fc_dim=fc_dim):
        super(LSTM_Attn, self).__init__()
        self.lstm = nn.LSTM(8, lstm_dim, batch_first=True, bidirectional=True, num_layers=lstm_layers)
        self.attention = nn.Sequential(
                        nn.Linear(lstm_dim*2, da),
                        nn.ReLU(),
                        nn.Linear(da, r)
                        )
        self.softmax = nn.Softmax(dim=1)
        self.fc1 = nn.Linear(lstm_dim*2*r, fc_dim)
        self.fc2 = nn.Linear(fc_dim, out_dim)

    def forward(self, inputs, hidden0=None):
        out, (hidden_state, cell_state) = self.lstm(inputs, hidden0)
        attn = self.softmax(self.attention(out))
        m = [0] * r
        for a in range(r):
            m[a] = (out*attn[:,:,a].unsqueeze(2)).sum(dim=1)
        out = torch.cat(m, dim=1)
        out = self.fc1(out)
        out = self.fc2(out)
        # out = self.softmax(out)
        return out, attn


# class LSTM(nn.Module):
#     def __init__(self, lstm_dim=lstm_dim, lstm_layers=lstm_layers, out_dim=2):
#         super(LSTM, self).__init__()
#         self.lstm = nn.LSTM(8, lstm_dim, batch_first=True, bidirectional=True, num_layers=lstm_layers)
#         self.fc1 = nn.Linear(lstm_dim*2, 128)
#         self.fc2 = nn.Linear(128, out_dim)
#         self.relu = nn.ReLU()
#         self.softmax = nn.Softmax(dim=1)

#     def forward(self, inputs, hidden0=None):
#         out, (hidden_state, cell_state) = self.lstm(inputs, hidden0)
#         out = torch.cat([hidden_state[0], hidden_state[1]], dim=1)
#         out = self.fc1(out)
#         out = self.fc2(out)
#         # out = self.softmax(out)
#         return out


# class CNN(nn.Module):
#     def __init__(self, dim=dim, channels=50, out_dim=2):
#         super(CNN, self).__init__()
#         self.conv1 = nn.Conv1d(in_channels=8, out_channels=channels, kernel_size=9, stride=1)
#         self.fc1 = nn.Linear((dim-9+1)*channels, 100)
#         self.drop = nn.Dropout(0.1)
#         self.fc2 = nn.Linear(100, out_dim)
#         self.relu = nn.ReLU()
#         self.softmax = nn.Softmax(dim=1)

#     def forward(self, inputs, hidden0=None):
#         out = torch.reshape(inputs, (inputs.shape[0], inputs.shape[2], inputs.shape[1]))
#         out = self.conv1(out)
#         out = out.view(out.size(0),-1)
#         out = self.fc1(out)
#         out = self.relu(out)
#         out = self.drop(out)
#         out = self.fc2(out)
#         # out = self.softmax(out)
#         return out


def GenomePred(genome_file, model_path, outfile, dim, batch_size, device):
    soft = nn.Softmax(dim=1)
    model = LSTM_Attn()
    model.load_state_dict(torch.load(model_path, map_location=device))
    model = model.to(device)
    ng_chrs = list('BDEFHIJKLMNOPQRSUVWXYZ0123456789-')
    with open(outfile, 'w') as fo:
        for genome in SeqIO.parse(genome_file, 'fasta'):
            node = genome.id
            seq = re.sub('[RMWSYKHBDV]', 'N', str(genome.seq).upper())
            for j in range(win+1, int(len(seq)-win), batch_size):
                l_f_num = []
                l_r_num = []
                l_posi = []
                for k in range(batch_size):
                    nucl = seq[j+k-win-1:j+k+win]
                    nucl = nucl.upper()
                    if (not 'N' in nucl) and (len(nucl) == dim):
                    # if not any((s in nucl) for s in ng_chrs) and (len(nucl) == dim):
                        f_num = onehot_plus(nucl, 'forward')
                        r_num = onehot_plus(nucl, 'revcom')
                        l_f_num.append(f_num)
                        l_r_num.append(r_num)
                        l_posi.append(j+k)

                if len(l_posi) > 0:
                    f_data = torch.tensor(l_f_num)
                    f_data = f_data.to(device)
                    f_output = soft(model(f_data)[0])   # Attn
                    # f_output = soft(model(f_data))  # LSTM & CNN
                    r_data = torch.tensor(l_r_num)
                    r_data = r_data.to(device)
                    r_output = soft(model(r_data)[0])   # Attn
                    # r_output = soft(model(r_data))  # LSTM & CNN

                    for l in range(len(l_posi)):
                        if f_output[l][0] > 0.5 or r_output[l][0] > 0.5:
                            # print('{}\t{}\t{:.4g}\t{:.4g}\n'.format(node, l_posi[l], f_output[l][0], r_output[l][0]), end='')
                            fo.write('{}\t{}\t{:.4g}\t{:.4g}\n'.format(node, l_posi[l], f_output[l][0], r_output[l][0]))
                            
                    del f_data, f_output, r_data, r_output
                    gc.collect()


GenomePred(genome_file=genome_file, model_path=model_path, outfile=outfile, dim=dim, batch_size=batch_size, device=device)

