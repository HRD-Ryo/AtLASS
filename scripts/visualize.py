#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
import sys
import numpy as np
import random
import torch
import torch.nn as nn


output = './model.onnx'


dim = 201
mer = 1
lstm_dim = 1024
lstm_layers = 1
da = 1024
fc_dim = 1024
batch_size = 200
lr = 0.0001
num_stop = 50


class LSTM_3Attn(nn.Module):
    def __init__(self, mer=mer, lstm_dim=lstm_dim, lstm_layers=lstm_layers, out_dim=2, da=da, r=3, fc_dim=fc_dim):
        super(LSTM_3Attn, self).__init__()
        self.lstm = nn.LSTM(4**mer+4, lstm_dim, batch_first=True, bidirectional=True, num_layers=lstm_layers)
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
        m1 = (out*attn[:,:,0].unsqueeze(2)).sum(dim=1)
        m2 = (out*attn[:,:,1].unsqueeze(2)).sum(dim=1)
        m3 = (out*attn[:,:,2].unsqueeze(2)).sum(dim=1)
        out = torch.cat([m1, m2, m3], dim=1)
        out = self.fc1(out)
        out = self.fc2(out)
        # out = self.softmax(out)
        return out, attn



model = LSTM_3Attn()
model.eval()

input = torch.randn(200, 201, 8)
torch.onnx.export(model, input, output)