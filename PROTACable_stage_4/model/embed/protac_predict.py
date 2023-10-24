"""
------------------------------------------------------------------------------------------
Copyright (C) 2023 Hazem Mslati. All Rights Reserved.

This script is provided "AS IS" without any warranty.
Unauthorized copying, modification, or distribution is prohibited.
------------------------------------------------------------------------------------------
"""

import os,re
import time
import argparse
import pytorch_lightning as pl
import joblib
from tqdm import tqdm
import pandas as pd

# import torch.multiprocessing as mp
# mp.set_start_method('spawn')

import torch
import torch.nn as nn
import torch.nn.functional as F

from pytorch_lightning.loggers import TensorBoardLogger
from pytorch_lightning.callbacks import ModelCheckpoint
from torch.utils import data
from torch.utils.data import DataLoader

from torch.nn.functional import sigmoid

import dgl
import numpy as np

from se3_transformer_protacs_nolstm_sep import *
import warnings
warnings.filterwarnings("ignore")
# try: torch.multiprocessing.set_start_method('spawn')
# except: pass
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


def _collate_fn(batch):
    nodes = []
    pairs = []
    bonds = []
    init_xyz = []
    init_pos = []
    init_atom = []
    init_CA = []
    label_xyz = []
    seq_l = []
    pdbs = []
    #res_atom = []
    for iter,item in enumerate(batch):
        nodes.append(item[0])
        init_xyz.append(item[1])
        init_pos.append(item[2])
        init_atom = item[3]
        init_CA.append(item[4])
        label_xyz.append(item[5])
        pdbs.append(item[6])
        pairs.append(item[7])
        bonds.append(item[8])
        l = item[1].shape[0]
        featurized = item[9] # [0] if item[9][0].shape else item[9] #.squeeze()
        sep = np.array(item[-1])
        seq_l.append(l)
    bsz = len(init_xyz)
    nodes = [torch.from_numpy(item).float() for item in nodes]
    nodes = torch.cat(nodes)
    pairs = [torch.from_numpy(item).float() for item in pairs]
    pairs = torch.cat(pairs)
    bonds = [torch.from_numpy(item).float() for item in bonds]
    bonds = torch.cat(bonds)
    init_xyz = [torch.from_numpy(item).float() for item in init_xyz]
    init_xyz = torch.cat(init_xyz)
    init_pos = [torch.from_numpy(item).float() for item in init_pos]
    init_pos = torch.cat(init_pos)
    init_CA = [torch.from_numpy(item).float() for item in init_CA]
    init_CA = torch.cat(init_CA)
    label_xyz = [torch.tensor(item) for item in label_xyz]
    label_xyz = torch.stack(label_xyz)
    features = [torch.tensor(item).float() for item in featurized]
    features = torch.stack(features)
    sep = [torch.from_numpy(sep)]
    sep = torch.stack(sep)


    return nodes, pairs, bonds, init_xyz, init_pos, init_atom, init_CA, label_xyz, seq_l, pdbs, features, sep


class MyDataset(data.Dataset):
    def __init__(self, file_list):
        self.file_list = file_list

    def __len__(self):
        return len(self.file_list)

    def __getitem__(self, idx):
        file = self.file_list[idx]
        try:
            with open(file, 'rb') as f:
                data = joblib.load(f)
            return data
        except EOFError:
            print(f"Error loading file: {file}")
            return None  # Or some other default value

def get_loader_from_file_list(file_list_path, num_workers, batch_size, shuffle=False):
    with open(file_list_path, 'r') as f:
        file_list = [line.strip() for line in f.readlines()]

    dataset = MyDataset(file_list)

    return DataLoader(dataset, shuffle=shuffle, pin_memory=True, num_workers=num_workers, batch_size=batch_size, collate_fn=_collate_fn)

def run_script(args):
    network = args.network
    test_seed = args.test_seed
    test_file = args.test_file
    num_workers = args.num_workers
    batch_size = args.batch_size
    checkpoint_path = args.checkpoint_path

    # Check that a valid checkpoint has been provided
    if checkpoint_path is None or not os.path.isfile(checkpoint_path):
        raise ValueError(f"Invalid checkpoint path: {checkpoint_path}")

    # Set the seed for reproducibility
    pl.seed_everything(seed=test_seed)

    test_loader = get_loader_from_file_list(test_file, num_workers, batch_size=1, shuffle=False)

    # Initialize the model and load from checkpoint
    test_model = globals()[network]()
    try: model = test_model.load_from_checkpoint(checkpoint_path)
    except: model = test_model.load_from_checkpoint(checkpoint_path, map_location=torch.device('cpu'))
    #checkpoint = torch.load(checkpoint_path)
    #model.load_state_dict(checkpoint['state_dict']) 

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model.to(device)

    # eval mode
    model.eval()

    val_step_outputs=[]
    val_step_targets=[]
    identifiers = []

    with torch.no_grad():
        for batch in tqdm(test_loader, desc="Processing batches"):

            nodes, pair, bond, init_pos, init_xyz, init_atom, init_CA, label, seq_l, pdbs, features, sep = batch
            identifiers.extend(pdbs)
            bsz = 1
            L = seq_l[0]
            print("processing batch ...")
            init_xyz = init_xyz.reshape(L, 1, 3)
            init_xyz = init_xyz.reshape(bsz,L,1,3)
            init_xyz = init_xyz.to(device)

            init_pos = init_pos.reshape(L, 1, 3)
            init_pos = init_pos.reshape(bsz,L,1,3)

            pair = pair.reshape(bsz,L,L,12)

            idx = torch.arange(L).long().view(1, L)
            idx = idx.to(device)
            idx = idx.reshape(bsz,L)
            idx = idx.to(device)
            
            nodes = model.norm_node1(nodes.unsqueeze(1))
            nodes = nodes.reshape(bsz,L,120)
            nodes = nodes.to(device)
            nodes = model.norm_node2(model.embed_node(nodes))

            pair = pair.permute(0,3,1,2)
            pair = model.norm_edge1(pair)
            pair = pair.permute(0,2,3,1)
            pair = pair.to(device)
            pair = model.norm_edge2(model.embed_e1(pair))

            rbf_feat = rbf(torch.cdist(init_xyz.reshape(bsz, L, 3), init_xyz.reshape(bsz, L, 3)))
            rbf_feat = rbf_feat.to(device)

            #pair = pair.to(device)

            bond = bond.reshape(1,L,L,1)
            bond = bond.to(device)

            pair = torch.cat((pair, rbf_feat, bond), dim=-1)
            pair = model.norm_edge3(model.embed_e2(pair))

            #nodes = nodes.to(device)
            pair = pair.to(device)
            #bond = bond.to(device)
            #init_pos = init_pos.to(device)
            #init_xyz = init_xyz.to(device)
            sep = sep.to(device)
            idx = idx.to(device)
            label = label.to(device)

            # define graph
            G = make_graph(init_xyz, pair, idx, sep, top_k=30)
            l1_feats = init_pos
            l1_feats = l1_feats.reshape(bsz*L,-1, 3)
            l1_feats = l1_feats.to(device)
            del idx, rbf_feat, pair, bond, init_pos, init_xyz, init_atom, init_CA, seq_l, pdbs

            # SE(3) Transformer
            shift = model.se3(G, nodes.reshape(bsz*L, -1, 1), l1_feats)
            output = shift['1'].reshape(bsz, L, -1)
            features = features.to(output.device)
            features = features.view(1, -1)
    
            output = model.AvgAttnConcatPoolingPredictor(output, features)

            # Make sure label_xyz is of shape (bsz,) and type LongTensor
            label = label.view(bsz).long()

            # Save outputs and targets per batch for later metric computations
            val_step_outputs.append(output.detach())
            val_step_targets.append(label.detach())

            del shift, output, nodes, G
            # TODO: other analysis?

    test_step_outputs_np = torch.sigmoid(torch.tensor(val_step_outputs)).numpy()

    if args.test_mode:
        # Remove the '_clean_interface_lframe.pdb' suffix from each identifier
        cleaned_identifiers = [id_str.replace('_clean_interface_lframe.pdb', '') for id_str in identifiers]
    
        df = pd.DataFrame({'Identifier': cleaned_identifiers, 'Predictions': test_step_outputs_np})
        df.to_csv('predictions_and_targets.csv', index=False)

    else: # to disk
        test_step_targets_np = np.array(val_step_targets)
        np.save('test_step_outputs_renewed_cont_sep_last.npy', test_step_outputs_np)
        np.save('test_step_targets_renewed_cont_sep_last.npy', test_step_targets_np)

def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('--network', type=str, required=False, default="SE3Network")
    parser.add_argument('--num_gpus', type=int, required=False, default=1)
    parser.add_argument('--num_workers', type=int, required=False, default=8)
    parser.add_argument('--batch_size', type=int, required=False, default=1)
    parser.add_argument('--test_seed', type=int, required=False, default=42)
    parser.add_argument('--test_mode', action='store_true', default=True)
    parser.add_argument('--test_file', type=str, required=True)
    parser.add_argument('--checkpoint_path', type=str, required=True, default=None)
#    parser.add_argument('--model_path', type=str, required=True, default=None)
    parser.add_argument('--device', type=str, required=False, default='cpu')


    args = parser.parse_args()

    return args

if __name__ == '__main__':
    args = parse_arguments()
    run_script(args)

