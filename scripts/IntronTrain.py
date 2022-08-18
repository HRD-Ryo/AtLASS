#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"    #### DDBJ
os.environ["CUDA_VISIBLE_DEVICES"]="0"    #### DDBJ
os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'max_split_size_mb:128'
import subprocess
import sys
from Bio import SeqIO
import numpy as np
import gc
import random
import torch
import torch.nn as nn
from sklearn.metrics import confusion_matrix
import re


####    IntronTrain with positive negative tsv
####    positive date = positive
####    negative data = positive revcom + negative


genome_file = sys.argv[1]
posi_file = sys.argv[2]
nega_file = sys.argv[3]
model_dir = sys.argv[4]
if sys.argv[5] in ['True', 'TRUE', 'true']:
    gpu = True
else:
    gpu = False

device = torch.device('cuda' if torch.cuda.is_available() and gpu else 'cpu')
print('# {} is available'.format(device))

dim = 201
win = int((dim-1)/2)
num_stop = 30
####    fast mode
if sys.argv[6] in ['True', 'TRUE', 'true']:
    lstm_dim = 256
    lstm_layers = 1
    da = 128
    r = 1
    fc_dim = 32
    lr = 0.0001
    batch_size = 1000
else:
    lstm_dim = 1024
    lstm_layers = 1
    da = 1024
    r = 3
    fc_dim = 1024
    lr = 0.0001
    batch_size = 500

print('# dim={}, lstm={}*{}, r={}, da={}, fc_dim={}, batch_size={}, lr={}'.format(dim, lstm_dim, lstm_layers, r, da, fc_dim, batch_size, lr))


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


def DataRead(in1=genome_file, in2=posi_file, in3=nega_file, dim=dim):
    d_genome = {}    ####    dictionary of genome id and seq
    for genome in SeqIO.parse(in1, 'fasta'):
        d_genome[genome.id] = re.sub('[RMWSYKHBDV]', 'N', str(genome.seq).upper())

    l_x_posi = []
    l_y_posi = []
    l_x_nega = []
    l_y_nega = []
    # ng_chrs = list('BDEFHIJKLMNOPQRSUVWXYZ0123456789-')

    ####    positive data
    with open(in2, 'r') as f:
        for i in f:
            l_line = i.rstrip().split('\t')
            nucl = d_genome[l_line[0]][int(l_line[1])-win-1:int(l_line[1])+win]
            if (not 'N' in nucl) and (len(nucl) == dim):
            # if (not any((s in nucl) for s in ng_chrs)) and (len(nucl) == dim):
                if l_line[2] == 'forward':
                    num = onehot_plus(nucl, 'forward')
                    l_x_posi.append(num)
                    l_y_posi.append([1, 0])    #### positive
                    num = onehot_plus(nucl, 'revcom')
                    l_x_nega.append(num)
                    l_y_nega.append([0, 1])    #### negative
                elif l_line[2] == 'revcom':
                    num = onehot_plus(nucl, 'revcom')
                    l_x_posi.append(num)
                    l_y_posi.append([1, 0])    #### positive
                    num = onehot_plus(nucl, 'forward')
                    l_x_nega.append(num)
                    l_y_nega.append([0, 1])    #### negative

    ####    negative data
    with open(in3, 'r') as f:
        for i in f:
            l_line = i.rstrip().split('\t')
            nucl = d_genome[l_line[0]][int(l_line[1])-win-1:int(l_line[1])+win]
            if (not 'N' in nucl) and (len(nucl) == dim):
            # if (not any((s in nucl) for s in ng_chrs)) and (len(nucl) == dim):
                num = onehot_plus(nucl, l_line[2])
                l_x_nega.append(num)
                l_y_nega.append([0, 1])    #### negative

    l_x = l_x_posi + l_x_nega
    l_y = l_y_posi + l_y_nega

    posi_num = len(l_x_posi)
    nega_num = len(l_x_nega)
    ratio = len(l_y_nega) / len(l_y_posi)
    print('# positive:{}, negative:{}, ratio:{}'.format(posi_num, nega_num, ratio))
    del l_x_posi, l_y_posi, l_x_nega, l_y_nega
    gc.collect()
    return l_x, l_y, ratio


def DataSplit(l_x, l_y, fold, fold_num=5):
    seed = random.randrange(100)
    random.seed(seed)
    random.shuffle(l_x)
    random.seed(seed)
    random.shuffle(l_y)
    random.seed()
    sep1 = round(len(l_x)/fold_num*fold)
    sep2 = round(len(l_x)/fold_num*(fold+1))
    l_test_x = l_x[sep1:sep2]
    l_test_y = l_y[sep1:sep2]
    l_train_x = l_x[:sep1] + l_x[sep2:]
    l_train_y = l_y[:sep1] + l_y[sep2:]
    print('# all_data:{}, train_data:{}, test_data:{}'.format(len(l_x), len(l_train_x), len(l_test_x)))
    return l_train_x, l_train_y, l_test_x, l_test_y


class DataSet:
    def __init__(self, l_x, l_y):
        self.x = torch.tensor(l_x)
        self.y = torch.tensor(l_y)
    def __len__(self):
        return self.x.shape[0]
    def __getitem__(self, i):
        return self.x[i], self.y[i]


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


def test(tensor):
    val, idx = tensor.max(axis=1)
    return idx


def modelmake(l_x, l_y, ratio, fold, epochs_num=200, batch_size=batch_size, lr=lr, num_stop=num_stop, device=device):
    print('# {} is used'.format(device))
    l_train_x, l_train_y, l_test_x, l_test_y = DataSplit(l_x=l_x, l_y=l_y, fold=fold)
    datasets = DataSet(l_train_x, l_train_y)
    dataloader = torch.utils.data.DataLoader(datasets, batch_size=batch_size, shuffle=True, drop_last=True)
    datasets_test = DataSet(l_test_x, l_test_y)
    dataloader_test = torch.utils.data.DataLoader(datasets_test, batch_size=batch_size, shuffle=False, drop_last=False)
    test_size = len(l_test_y)
    del l_train_x, l_train_y, l_test_x, l_test_y
    gc.collect()
    model = LSTM_Attn()
    model.to(device)
    weights = torch.tensor([1., 1.]).to(device)
    # weights = torch.tensor([(ratio-1)/1+1, 1.]).to(device)
    print("# weights: {}".format(weights))
    criterion = nn.CrossEntropyLoss(weight=weights)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    best_f1 = 0.0
    best_epoch = 1

    print('# epoch start')
    for epoch in range(epochs_num):
    ####    training
        running_loss = 0.0
        model.train()
        for data, label in dataloader:
            data = data.to(device)
            label = label.to(device)
            optimizer.zero_grad()
            output, attn = model(data)    #### LSTM attn
            # output = model(data)    #### LSTM & CNN
            loss = criterion(output, label.argmax(axis=1))
            loss.backward()
            optimizer.step()
            running_loss += loss.item() * data.size(0)

    ####    test
        run_out = []
        run_label = []
        running_loss_test = 0.0
        accuracy = 0
        model.eval()
        for data, label in dataloader_test:
            data = data.to(device)
            label = label.to(device)
            output, attn = model(data)    #### LSTM attn
            # output = model(data)    #### LSTM & CNN
            loss = criterion(output, label.argmax(axis=1))
            running_loss_test += loss.item() * data.size(0)

            output2 = test(output).tolist()
            label2 = test(label).tolist()
            run_out += output2
            run_label += label2

    ####    validation
            accuracy += torch.count_nonzero((test(output)-test(label)==0))
        accuracy_rate = accuracy / test_size
        test_loss = running_loss_test/len(datasets_test)

        tn, fp, fn, tp = confusion_matrix(run_label, run_out, labels=[1, 0]).flatten()
        acc = (tp+tn) / (tp+fn+fp+tn)
        if tp == 0:
            pre, rec, f1 = 0, 0, 0
        else:
            pre = tp / (tp+fp)
            rec = tp / (tp+fn)
            f1 = 2*pre*rec / (pre+rec)

        print('epoch:{}, train_loss:{:.5g}, test_loss:{:.5g}, TP:{}, FN:{}, FP:{}, TN:{}, acc:{:.4g}, pre:{:.4g}, rec:{:.4g}, f1:{:.4g}'\
            .format(epoch+1, running_loss/len(datasets), test_loss, tp, fn, fp, tn, acc, pre, rec, f1))

        if best_f1 <= f1:
            best_f1 = f1
            best_epoch = epoch
            model_state_path = '{}epoch{}_model_state.pth'.format(model_dir, epoch+1)
            torch.save(model.state_dict(), model_state_path)  ####    save model
        if epoch - best_epoch > num_stop:   ####    train finish !
            break

    cmd1 = 'cp {0}epoch{1}_model_state.pth {0}epoch{1}_best_model_state.pth'.format(model_dir, best_epoch+1)
    cmd2 = 'rm {}epoch*[0-9]_model_state.pth'.format(model_dir)
    os.system(cmd1)
    os.system(cmd2)



l_x, l_y, ratio = DataRead()


modelmake(l_x=l_x, l_y=l_y, ratio=ratio, fold=0)

# for fold in range(5):
#     print('# fold:{}'.format(fold))
#     modelmake(l_x=l_x, l_y=l_y, ratio=ratio, fold=fold)

