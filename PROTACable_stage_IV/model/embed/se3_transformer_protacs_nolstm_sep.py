"""
------------------------------------------------------------------------------------------
Copyright (C) 2023 Hazem Mslati. All Rights Reserved.

This script is provided "AS IS" without any warranty.
Unauthorized copying, modification, or distribution is prohibited.
------------------------------------------------------------------------------------------
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import DataLoader
import torchvision.transforms as transforms
import pytorch_lightning as pl
import dgl
import math
import copy
import numpy as np
from typing import List, Tuple
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torchmetrics import F1Score
from torchmetrics import Accuracy
from sklearn.metrics import roc_auc_score

from sklearn.metrics import r2_score
from contextlib import nullcontext

from typing import Dict

import dgl
import dgl.function as fn
from dgl.nn.pytorch.softmax import edge_softmax
from dgl.nn.pytorch.glob import AvgPooling, MaxPooling

from packaging import version

from equivariant_attention.modules import get_basis_and_r, GSE3Res, GNormBias
from equivariant_attention.modules import GConvSE3, GNormSE3
from equivariant_attention.fibers import Fiber


network='SE3Network'

# test_model = globals()[network]()

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

SE3_param = {
    "num_layers": 2,  # from 2 to 1
    "num_channels": 16,  # from 16 to 12
    "num_degrees": 2,  # Keep this the same, as it may have a significant impact on performance.
    "l0_in_features": 32,  # from 32 to 24
    "l0_out_features": 8,  # from 8 to 4
    "l1_in_features": 1,  # Keep this the same, assuming input features cannot be changed
    "l1_out_features": 1,  # Keep this the same, assuming output features cannot be changed
    "num_edge_features": 32,  # from 32 to 24
    "div": 2,  # Assuming this is fixed and related to the internal workings of the SE3 transformer
    "n_heads": 8  # from 4 to 2
}

tfn_param = {
    "num_layers": 2,  # from 2 to 1
    "num_channels": 16,  # from 16 to 12
    "num_degrees": 2,  # Keep this the same, as it may have a significant impact on performance.
    "l0_in_features": 32,  # from 32 to 24
    "l0_out_features": 8,  # from 8 to 4
    "l1_in_features": 1,  # Keep this the same, assuming input features cannot be changed
    "l1_out_features": 1,  # Keep this the same, assuming output features cannot be changed
    "num_edge_features": 32  # from 32 to 24
}




class SE3Transformer(nn.Module):
    """SE(3) equivariant GCN with attention"""
    def __init__(self, num_layers=3, num_channels=32, num_degrees=3, n_heads=8, div=4,
                 si_m='1x1', si_e='att',
                 l0_in_features=32, l0_out_features=32,
                 l1_in_features=3, l1_out_features=3,
                 num_edge_features=32, x_ij=None):
        super().__init__()
        # Build the network
        self.num_layers = num_layers
        self.num_channels = num_channels
        self.num_degrees = num_degrees
        self.edge_dim = num_edge_features
        self.div = div
        self.n_heads = n_heads
        self.si_m, self.si_e = si_m, si_e
        self.x_ij = x_ij

        if l1_out_features > 0:
            fibers = {'in': Fiber(dictionary={0: l0_in_features, 1: l1_in_features}),
                           'mid': Fiber(self.num_degrees, self.num_channels),
                           'out': Fiber(dictionary={0: l0_out_features, 1: l1_out_features})}
        else:
            fibers = {'in': Fiber(dictionary={0: l0_in_features, 1: l1_in_features}),
                           'mid': Fiber(self.num_degrees, self.num_channels),
                           'out': Fiber(dictionary={0: l0_out_features})}

        blocks = self._build_gcn(fibers)
        self.Gblock = blocks

    def _build_gcn(self, fibers):
        # Equivariant layers
        Gblock = []
        fin = fibers['in']
        for i in range(self.num_layers):
            Gblock.append(GSE3Res(fin, fibers['mid'], edge_dim=self.edge_dim,
                                  div=self.div, n_heads=self.n_heads,
                                  learnable_skip=True, skip='cat',
                                  selfint=self.si_m, x_ij=self.x_ij))
            Gblock.append(GNormBias(fibers['mid']))
            fin = fibers['mid']
        Gblock.append(
            GSE3Res(fibers['mid'], fibers['out'], edge_dim=self.edge_dim,
                    div=1, n_heads=min(1, 2), learnable_skip=True,
                    skip='cat', selfint=self.si_e, x_ij=self.x_ij))
        return nn.ModuleList(Gblock)

    @torch.cuda.amp.autocast(enabled=False)
    def forward(self, G, type_0_features, type_1_features):
        # Compute equivariant weight basis from relative positions
        basis, r = get_basis_and_r(G, self.num_degrees-1)
        h = {'0': type_0_features, '1': type_1_features}
        for layer in self.Gblock:
            h = layer(h, G=G, r=r, basis=basis)
        return h


def rbf(D):
    # Distance radial basis function
    D_min, D_max, D_count = 0., 20., 36
    D_mu = torch.linspace(D_min, D_max, D_count).to(D.device)
    D_mu = D_mu[None,:]
    D_sigma = (D_max - D_min) / D_count
    D_expand = torch.unsqueeze(D, -1)
    RBF = torch.exp(-((D_expand - D_mu) / D_sigma)**2)
    return RBF


class LayerNorm(nn.Module):
    def __init__(self, d_model, eps=1e-5):
        super(LayerNorm, self).__init__()
        self.a_2 = nn.Parameter(torch.ones(d_model))
        self.b_2 = nn.Parameter(torch.zeros(d_model))
        self.eps = eps

    def forward(self, x):
        mean = x.mean(-1, keepdim=True)
        std = torch.sqrt(x.var(dim=-1, keepdim=True, unbiased=False) + self.eps)
        x = self.a_2*(x-mean)
        x /= std
        x += self.b_2
        return x


def make_graph(xyz, pair, idx, sep, top_k, kmin=12):

    B, L, A = xyz.shape[0], xyz.shape[1], xyz.shape[2]
    device = xyz.device

    top_k_value = top_k

    # Flatten residues and atoms into a single dimension
    xyz_flattened = xyz.view(B, L*A, 3)

    # Calculate pairwise distances for all atoms
    D = torch.cdist(xyz_flattened, xyz_flattened) + torch.eye(L*A, device=device).unsqueeze(0) * 999.9

    # Adjust sequence separation for all atoms
    idx_expanded = idx.unsqueeze(-1).repeat(1, 1, A).view(B, L*A)
    sep_adjusted = idx_expanded[:, None, :] - idx_expanded[:, :, None]
    sep_adjusted = sep_adjusted.abs() + torch.eye(L*A, device=device).unsqueeze(0) * 999.9

    # Get top_k neighbors
    D_neigh, E_idx = torch.topk(D, min(top_k_value, L*A), largest=False)
    topk_matrix = torch.zeros((B, L*A, L*A), device=device)
    topk_matrix.scatter_(2, E_idx, 1.0)

    cond = torch.logical_or(topk_matrix > 0.0, sep_adjusted < kmin)
    b, i, j = torch.where(cond)

    # Create a mask to identify protein-protein, ligand-ligand, and protein-ligand interactions
    protein_mask_i = i < sep
    protein_mask_j = j < sep
    ligand_mask_i = i >= sep
    ligand_mask_j = j >= sep

    protein_protein = protein_mask_i & protein_mask_j
    ligand_ligand = ligand_mask_i & ligand_mask_j
    protein_ligand = (protein_mask_i & ligand_mask_j) | (ligand_mask_i & protein_mask_j)

    # Assign weights based on the mask
    weights = torch.ones_like(i, dtype=torch.float32, device=device)
    weights[protein_protein] = 1.0
    weights[ligand_ligand] = 1.5
    weights[protein_ligand] = 3.0

    src = b * L*A + i
    tgt = b * L*A + j
    G = dgl.graph((src, tgt), num_nodes=B*L*A).to(device)
    G.edata['d'] = (xyz_flattened[b, j, :] - xyz_flattened[b, i, :]).detach()
    G.edata['w'] = pair[b, i // A, j // A] * weights.unsqueeze(-1)  # Adjust indexing if pair contains info for all atom combinations

    return G

class PairwiseConv(nn.Module):
    """SE(3)-equivariant convolution between two single-type features"""
    def __init__(self, degree_in: int, nc_in: int, degree_out: int,
                 nc_out: int, edge_dim: int=0):
        """SE(3)-equivariant convolution between a pair of feature types.

        This layer performs a convolution from nc_in features of type degree_in
        to nc_out features of type degree_out.

        Args:
            degree_in: degree of input fiber
            nc_in: number of channels on input
            degree_out: degree of out order
            nc_out: number of channels on output
            edge_dim: number of dimensions for edge embedding
        """
        super().__init__()
        # Log settings
        self.degree_in = degree_in
        self.degree_out = degree_out
        self.nc_in = nc_in
        self.nc_out = nc_out

        # Functions of the degree
        self.num_freq = 2*min(degree_in, degree_out) + 1
        self.d_out = 2*degree_out + 1
        self.edge_dim = edge_dim

        # Radial profile function
        self.rp = RadialFunc(self.num_freq, nc_in, nc_out, self.edge_dim)

    def forward(self, feat, basis):
        # Get radial weights
        R = self.rp(feat)
        kernel = torch.sum(R * basis[f'{self.degree_in},{self.degree_out}'], -1)
        return kernel.view(kernel.shape[0], self.d_out*self.nc_out, -1)

class AttentionPool1d(nn.Module):
    "Attention for Learned Aggregation"
    def __init__(self, ni, bias=True, norm=nn.LayerNorm):
        super().__init__()
        self.norm = norm(ni)
        self.q = nn.Linear(ni, ni, bias=bias)
        self.vk = nn.Linear(ni, ni*2, bias=bias)
        self.proj = nn.Linear(ni, ni)

    def forward(self, x, cls_q):
        x = self.norm(x)
        B, N, C = x.shape

        q = self.q(cls_q)
        k, v = self.vk(x).reshape(B, N, 2, C).permute(2, 0, 1, 3).chunk(2, 0)

        attn = (q @ k.transpose(-2, -1))
        attn = attn.softmax(dim=-1)

        x = (attn @ v).transpose(1, 2).reshape(B, C)
        return self.proj(x)

class AvgAttnConcatPooling1d(nn.Module):
    def __init__(self, ni, se3_dim, n_out, attn_bias=True, ffn_expand=3, norm=nn.LayerNorm, act_cls=nn.GELU):
        super().__init__()
        self.cls_q = nn.Parameter(torch.zeros([1,ni]))
        self.attn = AttentionPool1d(ni, attn_bias, norm)
        self.norm1 = norm(ni)
        self.norm2 = norm(ni)
        self.norm5 = norm(2248)
        self.ln = nn.LayerNorm(se3_dim + ni)
        self.ffn = nn.Sequential(
            nn.Linear(ni, int(ni*ffn_expand)),
            act_cls(),
            norm(int(ni*ffn_expand)),
            nn.Linear(int(ni*ffn_expand), ni)
        )
        self.norm3 = norm(ni)
        self.out_transform = TCN()  # OutTransform to reduce feature and sequence dimensions
        self.norm4 = norm(1024)  # New number of features after OutTransform
        self.act = act_cls()
        nn.init.trunc_normal_(self.cls_q, std=0.02)
        #self.apply(self._init_weights)
        # Fully connected layers for regression
        self.fc1 = nn.Linear(se3_dim + ni + 2248, 512)
        self.fc2 = nn.Linear(512, 128)
        self.fc3 = nn.Linear(128, 64)
        self.regression = nn.Linear(64, n_out)
        self.dropout = nn.Dropout(0.3)

    def forward(self, x, ext_feats):
        a = self.cls_q + self.norm1(self.attn(x, self.cls_q))
        a = a + self.ffn(self.norm2(a))
        x_pooled = self.out_transform(x)  # Apply OutTransform
        x = self.act(torch.cat([self.norm4(x_pooled), self.norm3(a), self.norm5(ext_feats)], dim=1))
        x = self.dropout(x)  # Apply dropout
        x = self.act(self.fc1(x))  # Shape: [B, 512], replaced ReLU with GELU
        x = self.dropout(x)  # Apply dropout
        x = self.act(self.fc2(x))  # Shape: [B, 128], replaced ReLU with GELU
        x = self.dropout(x)
        x = self.act(self.fc3(x))  # Shape: [B, 64], replaced ReLU with GELU
        x = self.dropout(x)
        x = self.regression(x)  # Shape: [B, n_out]
        return x

  #  @torch.no_grad()
  #  def _init_weights(self, m):
  #      if isinstance(m, nn.Linear):
  #          nn.init.trunc_normal_(m.weight, std=0.02)
  #          if m.bias is not None:
  #              nn.init.constant_(m.bias, 0)


class TCNBlock(nn.Module):
    def __init__(self, in_channels, out_channels, kernel_size=3, dropout=0.3):
        super(TCNBlock, self).__init__()
        self.conv = nn.Conv1d(in_channels, out_channels, kernel_size, dilation=2)
        self.dropout = nn.Dropout(dropout)
        self.act = nn.GELU()
        self.downsample = nn.Conv1d(in_channels, out_channels, 1) if in_channels != out_channels else None

    def forward(self, x):
        # x: (batch_size, in_channels, L)
        residual = x
        if self.downsample is not None:
            residual = self.downsample(residual)  # (batch_size, out_channels, L)

        padding = (self.conv.weight.size(2) - 1) * self.conv.dilation[0]

        x = F.pad(x, (padding, 0))  # (batch_size, in_channels, L + padding)

        x = self.act(self.conv(x))  # (batch_size, out_channels, L)
        x = self.dropout(x)  # (batch_size, out_channels, L)
        x += residual  # (batch_size, out_channels, L)
        return x

class TCN(nn.Module):
    def __init__(self):
        super(TCN, self).__init__()
        self.tcn1 = TCNBlock(3, 64)
        self.tcn2 = TCNBlock(64, 32)
        self.tcn3 = TCNBlock(32, 1)

        self.adaptive_pool = nn.AdaptiveAvgPool1d(1024)

#        self.lstm = nn.LSTM(input_size=1, hidden_size=512, num_layers=2, batch_first=True, bidirectional=True, dropout=0.3)

    def forward(self, x):
        # x: (batch_size, L, 3)
        x = x.permute(0, 2, 1)  # x: (batch_size, 3, L)
        x = self.tcn1(x)  # x: (batch_size, 64, L)
        x = self.tcn2(x)  # x: (batch_size, 32, L)
        x = self.tcn3(x)  # x: (batch_size, 1, L)

        # apply adaptive pooling on TCN output
        pooled_x = self.adaptive_pool(x).view(x.size(0), -1)  # pooled_x: (batch_size, 1024)

 #       x = x.permute(0, 2, 1)  # x: (batch_size, L, 1)
  #      output, (hidden_state, cell_state) = self.lstm(x)  # x: (batch_size, L, 2*hidden_size); hidden_state: (2*num_layers, batch_size, hidden_size)

   #     lstm_output = torch.cat((hidden_state[-2,:,:], hidden_state[-1,:,:]), dim=1)  # lstm_output: (batch_size, 2*hidden_size)

        # concat pooled TCN output and LSTM output
        #x = torch.cat((pooled_x, lstm_output), dim=1)  # x: (batch_size, 1024 + 2*hidden_size)
        return pooled_x #, lstm_output


class MultiHeadAttention(nn.Module):
    def __init__(self, num_heads=8, d_model=64, dropout=0.1):
        super().__init__()
        assert d_model % num_heads == 0, "num_heads must be a divisor of d_model"
        
        self.d_model = d_model
        self.num_heads = num_heads
        self.head_dim = d_model // num_heads
        
        self.qkv = nn.Linear(d_model, d_model * 3, bias=False)  # queries, keys, values
        self.o = nn.Linear(d_model, d_model, bias=False)  # output linear layer
        self.dropout = nn.Dropout(dropout)

    def forward(self, x, mask=None):
        batch_size, seq_len, _ = x.size()
        
        qkv = self.qkv(x).view(batch_size, seq_len, 3, self.num_heads, self.head_dim).permute(2, 0, 3, 1, 4)
        q, k, v = qkv[0], qkv[1], qkv[2]  # [batch_size, num_heads, seq_len, head_dim]
        
        attn_scores = torch.matmul(q, k.transpose(-2, -1)) / math.sqrt(self.head_dim + 1e-9)
        
        if mask is not None:
            attn_scores = attn_scores.masked_fill(mask == 0, -1e9)
            
        attn_probs = F.softmax(attn_scores, dim=-1)
        attn_probs = self.dropout(attn_probs)
        
        attn_output = torch.matmul(attn_probs, v)  # [batch_size, num_heads, seq_len, head_dim]
        attn_output = attn_output.transpose(1, 2).contiguous().view(batch_size, seq_len, -1)  # [batch_size, seq_len, d_model]
        
        return self.o(attn_output)  # [batch_size, seq_len, d_model]
    
class AttentionPooling(nn.Module):
    def __init__(self, input_dim, attention_dim, output_dim):
        super().__init__()
        self.attention_proj = nn.Linear(input_dim, attention_dim)
        self.relu = nn.ReLU()
        self.attention_softmax = nn.Softmax(dim=1)
        self.fc = nn.Linear(input_dim, output_dim)

    def forward(self, x):
        # Compute attention scores
        attn_scores = self.attention_proj(x)
        attn_scores = self.relu(attn_scores)
        attn_scores = self.attention_softmax(attn_scores)

        # Multiply features by attention scores and sum to get pooled features
        pooled = torch.sum(attn_scores * x, dim=1)

        out = self.fc(pooled)
        
        return out

class RMSE(nn.Module):
    def __init__(self):
        super().__init__()
        self.mse = nn.MSELoss()

    def forward(self, yhat, y):
        return torch.sqrt(self.mse(yhat, y))

class PositionalEncoding(nn.Module):
    def __init__(self, d_model, dropout=0.1, max_len=5000):
        super(PositionalEncoding, self).__init__()
        self.dropout = nn.Dropout(p=dropout)
        self.d_model = d_model
        self.padding = None
        if d_model % 2 != 0:  # add padding if d_model is odd
            self.d_model += 1
            self.padding = nn.ConstantPad1d((0, 1), 0)

        # Compute the positional encodings once in log space.
        pe = torch.zeros(max_len, self.d_model)
        position = torch.arange(0, max_len).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, self.d_model, 2) *
                             -(math.log(10000.0) / self.d_model))
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0)
        self.register_buffer('pe', pe)

    def forward(self, x):
        if self.padding is not None:
            x = self.padding(x)
        x = x + self.pe[:, :x.size(1)]
        if self.padding is not None:
            x = x[..., :self.d_model-1]  # remove the dummy dimension
        return self.dropout(x)


class Autoencoder(nn.Module):
    def __init__(self, input_dim, embedding_dim):
        super(Autoencoder, self).__init__()
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, input_dim//2),
            nn.ReLU(),
            nn.Linear(input_dim//2, embedding_dim),
        )

    def forward(self, x):
        x = self.encoder(x)
        return x

class SE3Network(pl.LightningModule):
    def __init__(self, *args, **kwargs):
        super(SE3Network, self).__init__()
        self.learning_rate = 0.0005

        #self.se3 = TFN(**tfn_param)
        self.se3 = SE3Transformer(**SE3_param)
#         self.lossfn = nn.MSELoss(reduction='sum')
        self.norm_node1 = nn.InstanceNorm1d(120)
        self.embed_node = nn.Linear(120, SE3_param['l0_in_features'])
        self.norm_node2 = LayerNorm(SE3_param['l0_in_features'])
        self.norm_edge1 = nn.InstanceNorm2d(12)
        self.embed_e1 =  nn.Linear(12, SE3_param['num_edge_features'])
        self.embed_e2 = nn.Linear(SE3_param['num_edge_features']+37, SE3_param['num_edge_features'])
        self.norm_edge2 = LayerNorm(SE3_param['num_edge_features'])
        self.norm_edge3 = LayerNorm(SE3_param['num_edge_features'])

        # Attributes to save batch outputs and targets
        self.training_step_outputs = []
        self.training_step_targets = []
        self.val_step_outputs = []
        self.val_step_targets = []
        self.val_loss=[]

        self.f1 = F1Score(num_classes=2, average='macro', task='binary')
        self.accuracy = Accuracy(task="binary")
        
        # Set d_model to match the dimension of the data you are working with
        # Set d_model to match the dimension of the data you are working with
        self.d_model = 1024 
        #self.autoencoder_shift_output = Autoencoder(3, self.d_model // 2)  # Apply autoencoder to SE3 output
        #self.pe = PositionalEncoding(3, dropout=0.1)
        #self.attention = MultiHeadAttention(num_heads=1, d_model=3)
        #self.fc_module = FCModule(3, self.d_model)
        #self.attention_pooling = AttentionPooling(input_dim=3, attention_dim=3, output_dim=64)
        #self.features_autoencoder = Autoencoder(2248, self.d_model * 2)  # Apply autoencoder to features
        #self.TransformOutput = OutTransform()
        self.AvgAttnConcatPoolingPredictor = AvgAttnConcatPooling1d(ni=3, se3_dim=1024, n_out=1)

        self.mse = nn.MSELoss()
        self.lossfn = nn.BCEWithLogitsLoss()
        self.rmse = RMSE()


    def training_step(self, batch, batch_idx):
        #free, total = cuda.mem_get_info()
        #print("Memory free loop begin: ", free/(1024**2), "MB")
        nodes, pair, bond, init_pos, init_xyz, init_atom, init_CA, label, seq_l, pdbs, features, sep = batch
        bsz = 1
        L = seq_l[0]
        print("batch L: ", L)

        init_xyz = init_xyz.reshape(L, 1, 3)
        init_xyz = init_xyz.reshape(bsz,L,1,3)

        init_pos = init_pos.reshape(L, 1, 3)
        init_pos = init_pos.reshape(bsz,L,1,3)

        pair = pair.reshape(bsz,L,L,12)

        idx = torch.arange(L).long().view(1, L)
        idx = idx.to(device)
        idx = idx.reshape(bsz,L)

        nodes = self.norm_node1(nodes.unsqueeze(1))
        nodes = nodes.reshape(bsz,L,120)
        nodes = self.norm_node2(self.embed_node(nodes))

        pair = pair.permute(0,3,1,2)
        pair = self.norm_edge1(pair)
        pair = pair.permute(0,2,3,1)
        pair = self.norm_edge2(self.embed_e1(pair))

        rbf_feat = rbf(torch.cdist(init_xyz.reshape(bsz, L, 3), init_xyz.reshape(bsz, L, 3)))
        #rbf_feat = rbf_feat.to(device)

        # Moved tensors to device
        #pair = pair.to(device)
        
        bond = bond.reshape(1,L,L,1)
        #bond = bond.to(device)
        
        pair = torch.cat((pair, rbf_feat, bond), dim=-1)
        
        #nodes = nodes.to(device)
        #init_pos = init_pos.to(device)
        #init_xyz = init_xyz.to(device)
        #label = label.to(device)

        pair = self.norm_edge3(self.embed_e2(pair))

        # define graph
        G = make_graph(init_xyz, pair, idx, sep, top_k=50) 
        l1_feats = init_pos
        l1_feats = l1_feats.reshape(bsz*L,-1, 3)

        del idx, rbf_feat, pair, bond, init_pos, init_xyz, init_atom, init_CA, seq_l
        # SE(3) Transformer
       
        shift = self.se3(G, nodes.reshape(bsz*L, -1, 1), l1_feats)
        output = shift['1'].reshape(bsz, L, -1)

        #output = self.pe(output) 
        #output = self.attention(output)
        #output = self.TransformOutput(output)
        #output = self.attention_pooling(output)
        #output = output.squeeze(-1)
        #print("output : ", output.detach().cpu().numpy().shape)
        #mask = torch.isnan(features)
        #if torch.any(mask):
        #    with open('nan_pdbs.txt', 'a') as f:
        #        f.write(f'NaN detected in features. Corresponding pdb: {pdbs}\n')
        #features = torch.where(mask, torch.zeros_like(features), features)
        #features = features.to(output.device)

        # Apply autoencoder to features
        #features_transformed = self.features_autoencoder(features.squeeze(1))

        # Concatenate features to the output
        #concatenated_output = torch.cat([output, features], dim=-1)
        #concatenated_output = concatenated_output.unsqueeze(1)
        #print("concat : ", concatenated_output.detach().cpu().numpy().shape)

#        print("output_concat_attn: ", output.detach().cpu().numpy().shape)

        EPSILON = 1e-7
        features = features.to(output.device)
        features = features.view(1, -1)

        output = self.AvgAttnConcatPoolingPredictor(output, features)

        #probabilities = torch.nn.functional.softmax(output, dim=1)
        # Make sure label_xyz is of shape (bsz,) and type LongTensor
        label = label.view(bsz).long()

        #l1_regularization = torch.tensor(0., device=self.device)  # Use the same device as your parameters
        #for param in self.parameters():
        #    l1_regularization += torch.norm(param, 1)

        # Add L1 regularization to the loss
        loss = self.lossfn(output.view(-1), label.float())
        #reg_lambda = 1e-5
        #loss += reg_lambda * l1_regularization

        print('step train loss = ', loss.detach())

        #self.training_step_outputs.extend(y_pred)
        #self.training_step_targets.extend(y_true)

        batch_dictionary = {'loss': loss}
        self.log('train_loss', loss, on_step=False, on_epoch=True)
        del shift, output, loss, nodes, G
        #torch.cuda.empty_cache()
        #free, total = cuda.mem_get_info()
        #print("Memory free loop end: ", free/(1024**2), "MB")
        return batch_dictionary        

    
    def on_train_epoch_end(self):
        pass
        #train_all_outputs = self.training_step_outputs
        #train_all_targets = self.training_step_targets
        #f1_macro_epoch = f1_score(train_all_outputs, train_all_targets, average='macro')
        #self.log("training_f1_epoch", f1_macro_epoch, on_step=False, on_epoch=True)

        # free up the memory
        #self.training_step_outputs.clear()
        #self.training_step_targets.clear()
        

    def validation_step(self, batch, batch_idx):
        nodes, pair, bond, init_pos, init_xyz, init_atom, init_CA, label, seq_l, pdbs, features, sep = batch
        bsz = 1
        L = seq_l[0]

        init_xyz = init_xyz.reshape(L, 1, 3)
        init_xyz = init_xyz.reshape(bsz,L,1,3)

        init_pos = init_pos.reshape(L, 1, 3)
        init_pos = init_pos.reshape(bsz,L,1,3)

        pair = pair.reshape(bsz,L,L,12)

        idx = torch.arange(L).long().view(1, L)
        idx = idx.to(device)
        idx = idx.reshape(bsz,L)

        nodes = self.norm_node1(nodes.unsqueeze(1))
        nodes = nodes.reshape(bsz,L,120)
        nodes = self.norm_node2(self.embed_node(nodes))

        pair = pair.permute(0,3,1,2)
        pair = self.norm_edge1(pair)
        pair = pair.permute(0,2,3,1)
        pair = self.norm_edge2(self.embed_e1(pair))

        rbf_feat = rbf(torch.cdist(init_xyz.reshape(bsz, L, 3), init_xyz.reshape(bsz, L, 3)))
        rbf_feat = rbf_feat.to(device)

        pair = pair.to(device)

        bond = bond.reshape(1,L,L,1)
        bond = bond.to(device)

        pair = torch.cat((pair, rbf_feat, bond), dim=-1)
        pair = self.norm_edge3(self.embed_e2(pair))

        nodes = nodes.to(device)
        pair = pair.to(device)
        bond = bond.to(device)
        init_pos = init_pos.to(device)
        init_xyz = init_xyz.to(device)
        idx = idx.to(device)
        label = label.to(device)

        # define graph
        G = make_graph(init_xyz, pair, idx, sep, top_k=50)
        l1_feats = init_pos
        l1_feats = l1_feats.reshape(bsz*L,-1, 3)
        del idx, rbf_feat, pair, bond, init_pos, init_xyz, init_atom, init_CA, seq_l, pdbs

        # SE(3) Transformer
        shift = self.se3(G, nodes.reshape(bsz*L, -1, 1), l1_feats)
        output = shift['1'].reshape(bsz, L, -1)

        #output = self.pe(output)

        #output = self.attention(output)
        #output = self.TransformOutput(output)
        #output = self.attention_pooling(output)
        #output = output.squeeze(-1)

        #features = features.to(output.device)

        # Apply autoencoder to features
        #features_transformed = self.features_autoencoder(features.squeeze(1))

        # Concatenate features to the output
        #concatenated_output = torch.cat([output, features], dim=-1)
        #concatenated_output = concatenated_output.unsqueeze(1)
        features = features.to(output.device)
        features = features.view(1, -1)

        output = self.AvgAttnConcatPoolingPredictor(output, features)
        # Make sure label_xyz is of shape (bsz,) and type LongTensor
        label = label.view(bsz).long()

        loss = self.lossfn(output.view(-1), label.float())

        # Make sure label_xyz is of shape (bsz,) and type LongTensor
        label = label.view(bsz).long()

        # Save outputs and targets per batch for later metric computations
        self.val_step_outputs.append(output.detach())
        self.val_step_targets.append(label.detach())

        if len(loss.shape) == 0:  # loss is a scalar
            loss = loss.view(1)  # Converts loss to a 1D tensor
        self.val_loss.extend(loss)

        batch_dictionary = {'val loss': loss}
        self.log('val_loss', loss, on_step=False, on_epoch=True)
        del shift, output, loss, nodes, G
        return batch_dictionary

    def on_validation_epoch_end(self):
        print("on validation epoch end")
        # Concatenate all outputs and targets
        all_outputs = torch.cat(self.val_step_outputs, dim=0)
        all_targets = torch.cat(self.val_step_targets, dim=0)

        # Compute other metrics
        y_pred_t = torch.sigmoid(all_outputs) > 0.5 # for binary classification
        y_pred_t = y_pred_t.view(-1) # reshape to match the shape of all_targets

        # Convert tensor to numpy array for sklearn
        all_outputs_np = all_outputs.cpu().numpy()
        all_targets_np = all_targets.cpu().numpy()

        print("outputs= ", all_outputs_np)
        print("targets= ", all_targets_np)
        # Calculate AUC-ROC
        auc_roc = roc_auc_score(all_targets_np, all_outputs_np)
        self.log('val_auc_roc', auc_roc, on_step=False, on_epoch=True, prog_bar=True, logger=True)

        avg_acc = self.accuracy(y_pred_t, all_targets)
        self.log('val_acc', avg_acc, on_step=False, on_epoch=True, prog_bar=True, logger=True)

        avg_f1 = self.f1(y_pred_t, all_targets)
        self.log('val_f1', avg_f1, on_step=False, on_epoch=True, prog_bar=True, logger=True)

        avg_loss = torch.stack([x for x in self.val_loss]).mean()
        self.log('val_loss', avg_loss, on_step=False, on_epoch=True, prog_bar=True, logger=True)



        # free up the memory
        self.val_step_outputs.clear()
        self.val_step_targets.clear()
        self.val_loss.clear()

   # def configure_optimizers(self):
    #    optimizer = torch.optim.RMSprop(self.parameters(), lr=1e-2, weight_decay=0.001)
     #   return [optimizer]

    def configure_optimizers(self):
        optimizer = torch.optim.Adam(self.parameters(), lr=1e-4, betas=(0.9, 0.999), weight_decay=0.005)
        return [optimizer]
