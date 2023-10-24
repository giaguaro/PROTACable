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

# import torch.multiprocessing as mp
# mp.set_start_method('spawn')

import torch
import torch.nn as nn
import torch.nn.functional as F

from pytorch_lightning.loggers import TensorBoardLogger
from pytorch_lightning.callbacks import ModelCheckpoint
from torch.utils import data
from torch.utils.data import DataLoader


import dgl
import numpy as np

from se3_transformer_protacs_nolstm import *
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
        #res_atom = item[9]
        l = item[1].shape[0]
        featurized = item[9][0] if item[9][0].shape else item[9] #.squeeze()
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

    return nodes, pairs, bonds, init_xyz, init_pos, init_atom, init_CA, label_xyz, seq_l, pdbs, features #, res_atom


# Define a new DataLoader function

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
    # parameters
    network = args.network
    out_path = args.out_path
    num_gpus = args.num_gpus
    num_workers = args.num_workers
    epochs = args.epochs
    batch_size = args.batch_size
    test_seed = args.test_seed
    debug_mode = args.debug_mode
    test_mode = args.test_mode
    train_file = args.train_file
    val_file = args.val_file
    test_file = args.test_file

    # Get the checkpoint path. If it's None or the file does not exist, the training will start from scratch
    checkpoint_path = args.checkpoint_path if args.checkpoint_path and os.path.isfile(args.checkpoint_path) else None

    pl.seed_everything(seed=test_seed)

    train_loader = get_loader_from_file_list(train_file, num_workers, batch_size, shuffle=True)
    val_loader = get_loader_from_file_list(val_file, num_workers, batch_size, shuffle=False)
    test_loader = get_loader_from_file_list(test_file, num_workers, batch_size=1, shuffle=False)

    start_time = time.time()

    test_model = globals()[network]()
    for p in test_model.parameters():
        if p.dim() > 1:
            nn.init.kaiming_uniform_(p)


    checkpoint_callback = ModelCheckpoint(
        monitor='val_acc',
        dirpath=out_path,
        filename=network + '-{epoch:02d}-{val_acc:.3f}-{val_f1:.3f}-{val_loss:.3f}',
        save_top_k=20,
        mode='max',
        every_n_epochs=1,
        save_last=True
    )

    if debug_mode:
        if not test_mode:
            trainer = pl.Trainer(precision=16, max_epochs=epochs, check_val_every_n_epoch=1,
                                 accumulate_grad_batches=batch_size,
                                 default_root_dir=out_path,
                                 accelerator='cuda',
                                 callbacks=[checkpoint_callback],
                                 resume_from_checkpoint=checkpoint_path,
                                 )
    else:
        logger = TensorBoardLogger(out_path, name="log")
        trainer = pl.Trainer(max_epochs=epochs, precision=16, check_val_every_n_epoch=1,
                             accumulate_grad_batches=batch_size*8,
                             default_root_dir=out_path,
                             accelerator='cuda',
                             strategy='ddp_find_unused_parameters_true',
                             logger=logger,
                             callbacks=[checkpoint_callback],
                             num_sanity_val_steps=0,
                            # resume_from_checkpoint=checkpoint_path,
                             )

    time1 = time.time()
    trainer.fit(test_model, train_loader, val_loader, ckpt_path=checkpoint_path)
    time2 = time.time()
    print('{} epochs takes {} seconds using {} GPUs.'.format(epochs, time2 - time1, num_gpus))

def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('--network', type=str, required=True)
    parser.add_argument('--out_path', type=str, required=True)
    parser.add_argument('--num_gpus', type=int, required=True)
    parser.add_argument('--num_workers', type=int, required=True)
    parser.add_argument('--epochs', type=int, required=True)
    parser.add_argument('--batch_size', type=int, required=True)
    parser.add_argument('--test_seed', type=int, required=True)
    parser.add_argument('--debug_mode', action='store_true')
    parser.add_argument('--test_mode', action='store_true')
    parser.add_argument('--train_file', type=str, required=True)
    parser.add_argument('--val_file', type=str, required=True)
    parser.add_argument('--test_file', type=str, required=True)
    parser.add_argument('--checkpoint_path', type=str, default=None) 

    args = parser.parse_args()

    return args

if __name__ == '__main__':
    args = parse_arguments()
    run_script(args)

