"""
------------------------------------------------------------------------------------------
Copyright (C) 2023 Hazem Mslati. All Rights Reserved.

This script is provided "AS IS" without any warranty.
Unauthorized copying, modification, or distribution is prohibited.
------------------------------------------------------------------------------------------
"""

import argparse
import gzip
import os, glob
from joblib import dump, load
import pickle
import warnings
from collections import defaultdict

import pandas as pd

import numpy as np
from Bio import SeqIO
from Bio.PDB import PDBParser, DSSP
from Bio.PDB.Polypeptide import is_aa

from rdkit import Chem
from rdkit.Chem import Mol

from plif_utilities import (RESIDUE_TYPES, pdb_to_fragments, clean_pdb,
                         get_atom_emb, merge_pdb_and_hetatm, set_lframe, SUBSTITUTIONS, COMMON_IONS, HETS, DNA_RNA)
from fragment_mol import FragmentPDB
from get_interface_ import nearest_residues_around_ligand
from align_filter import *
from compute_extra_feats import *


COMMON_IONS_UPPER = [ion.upper() for ion in COMMON_IONS]

STOP_ATOMS = {"C", "O", "N"}  # Stop checking if encountered a carbon, oxygen, or nitrogen

def get_seqs_from_pdb(pdb_file, updated):
    seq = []

    with open(pdb_file, 'r') as f:
        lines = f.readlines()

    prev_res_num = None
    res_to_append = None
    for line in lines:
        if line.startswith('ATOM') or (updated and line.startswith('HETATM')):
            res_num = line[22:26].strip()  # Extracting residue number
            residue_name = line[17:20].strip()
            atom_element = line[76:78].strip().upper()  # Extracting the element from the line
            atom_name = line[12:16].strip()

            # Append to seq if residue number changes
            if prev_res_num and prev_res_num != res_num:
                if res_to_append:
                    seq.append(res_to_append)
                res_to_append = None

            if atom_element in COMMON_IONS_UPPER:
                res_to_append = "ION"
                print("Found ion: ", atom_element)
            elif line.startswith('HETATM'):
                if not res_to_append:  # If not already determined as ION or other
                    res_to_append = "LIG"
            else:
                res_to_append = residue_name

            prev_res_num = res_num

    if res_to_append:  # To handle the last residue
        seq.append(res_to_append)

    print("Total number of ions to embedded is: ", seq.count('ION'))
    return seq

def get_seq_onehot(seq):
    seq_onehot = np.zeros((len(seq), len(RESIDUE_TYPES)))
    total_p = 0
    total_l = 0
 #   total_n = 0
    for i, res in enumerate(seq):
        # Check if the residue is part of the substitutions and update accordingly
        if res in SUBSTITUTIONS:
            res = SUBSTITUTIONS[res]

        if res not in RESIDUE_TYPES:
            res = "XXX"
        seq_onehot[i, RESIDUE_TYPES.index(res)] = 1
        if res not in HETS:
            total_p += 1
            seq_onehot[i, RESIDUE_TYPES.index("PROTEIN")] = 1
        if res == "LIG":
            total_l += 1
            seq_onehot[i, RESIDUE_TYPES.index("LIGAND")] = 1
    #    if res in DNA_RNA:
   #         total_n += 1
  #          seq_onehot[i, RESIDUE_TYPES.index("NUCLEIC")] = 1

    print("The total number of protein residues are: ", total_p)
    print("The total number of ligand residues are: ", total_l)
   # print("The total number of nucleic residues are: ", total_n)
    return seq_onehot


def get_rPos(seq):
    seq_len= len(seq)
    r_pos = np.linspace(0, 1, num=seq_len).reshape(seq_len, -1)
    return r_pos


def get_seq_feature(pdb_file, updated=False):
    seq = get_seqs_from_pdb(pdb_file, updated)
    node_feat = {
        'onehot': get_seq_onehot(seq),
        'rPosition': get_rPos(seq),
    }
    return node_feat, len(seq)


def filter_mols(side_chains_mols):
    filtered_mols = {k: list(v) for k, v in side_chains_mols.items()}
    for chain_name, mols in filtered_mols.items():
        for i in range(len(mols)):
            mol = mols[i]
            if mol is not None:
                try:
                    if (mol.GetNumHeavyAtoms() < 3 or
                        mol.GetAtomWithIdx(0).GetSymbol() != 'N' or
                        mol.GetAtomWithIdx(1).GetSymbol() != 'C' or
                        mol.GetAtomWithIdx(2).GetSymbol() != 'C' or
                        mol.GetAtomWithIdx(1).GetDegree() <= 1):
                        mols[i] = None
                except:
                    mols[i] = None
    return filtered_mols


def extract_res_atom_from_pdb(pdb_file_path):
    res_atom = defaultdict(list)
    with open(pdb_file_path, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                res_name = line[17:20].strip()
                atom_name = line[12:16].strip()
                res_atom[res_name].append(atom_name)
    return res_atom

def get_struc_feat(pdb_file, true_pdb, seq_len, allow_exceptions=False, protacs_embedding=False, ligand_docked=False):
    name = pdb_file.split('/')[-1].split('.')[0]
    side_chains_pdbs, ligands_pdbs = pdb_to_fragments(pdb_file, ligand_docked)
    side_chains_mols = {chain_name: [Chem.MolFromPDBBlock(residue) for residue in residues]
                        for chain_name, residues in side_chains_pdbs.items()}
    filtered_side_chain_mols = filter_mols(side_chains_mols)
    ligands_mols = {ligand_name: Chem.MolFromPDBBlock(ligand_pdb) for ligand_name, ligand_pdb in ligands_pdbs.items()}
    ligands_mols_lst = [Chem.MolFromPDBBlock(ligand_pdb) for _ , ligand_pdb in ligands_pdbs.items()]
    
    if any(ligand is not None for ligand in ligands_mols_lst):
        print("Detected ligands")
        print("Total ligands found: ", len(ligands_mols_lst))
        #print([[atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetSymbol() != 'H'] for mol in ligands_mols_lst if mol])
        fragment_pdb = FragmentPDB(protacs_embedding)
        accepted_fragments, fragment_counts, tot_orig_mols = fragment_pdb.split_and_write(ligands_mols_lst, f"{name}_ligand_lframe")
        updated_ligand = f"{name}_ligand_lframe.pdb"
        original_hets = f"{name}_ligand_lframe_original_mols.pdb"
        #print(fragment_counts)
    else:
        accepted_fragments, fragment_counts = None, None
        updated_ligand = None
        original_hets = None
        tot_orig_mols = 0
    
    protein_atom_embs, protein_atom_xyz, ligand_atom_embs, ligand_atom_xyz, protein_atom_lst, \
    ligand_atom_lst, CA_lst, protein_res_atom, ligand_res_atom, protein_covalent_bond, ligand_covalent_bond, \
    total_atoms, exclude_res = get_atom_emb(filtered_side_chain_mols, ligands_mols, fragment_counts, pdb_file, [1, seq_len], allow_exceptions)
    
    new_pdb_file = merge_pdb_and_hetatm(pdb_file, updated_ligand, f"{pdb_file.split('/')[-1].split('.')[0]}_lframe.pdb", exclude_res)
    merge_pdb_and_hetatm(pdb_file, original_hets, f"{pdb_file.split('/')[-1].split('.')[0]}_lframe_original.pdb", exclude_res)

    res_atom = extract_res_atom_from_pdb(f"{pdb_file.split('/')[-1].split('.')[0]}_lframe_original.pdb")

    node_edge_feat_rec = {}
    node_edge_feat_lig = {}
    combined_rec_edge_feats = {}
    
    node_edge_feat_rec['atom_emb'] = {
        'embedding': protein_atom_embs,
        'atom_pos': protein_atom_xyz,
        'atom_lst': protein_atom_lst,
        'CA_lst': CA_lst,
        'res_atom': protein_res_atom,
        'covalent': protein_covalent_bond
    }
    node_edge_feat_lig['atom_emb'] = {
        'embedding': ligand_atom_embs,
        'atom_pos': ligand_atom_xyz,
        'atom_lst': ligand_atom_lst,
        'CA_lst': CA_lst,
        'res_atom': ligand_res_atom,
        'covalent': ligand_covalent_bond
    }
    protein_res_atom.update(ligand_res_atom)
    combined_res_atom = protein_res_atom
    if ligand_covalent_bond is not None:
        print("also building covalent matrix for ligand")
        combined_covalent_bond = np.zeros((total_atoms[0] + total_atoms[1], total_atoms[0] + total_atoms[1]))
        combined_covalent_bond[:total_atoms[0], :total_atoms[0]] = protein_covalent_bond
        combined_covalent_bond[total_atoms[0]:, total_atoms[0]:] = ligand_covalent_bond
    else:
        combined_covalent_bond = protein_covalent_bond
    combined_rec_edge_feats['atom_emb'] = {
        'embedding': protein_atom_embs,
        'atom_pos': protein_atom_xyz,
        'atom_lst': protein_atom_lst + ligand_atom_lst,
        'CA_lst': CA_lst,
        'res_atom': [res_atom,tot_orig_mols],
        'covalent': combined_covalent_bond
    }
    return new_pdb_file, node_edge_feat_rec, node_edge_feat_lig, combined_rec_edge_feats

def get_item(pdb_file, true_pdb=None,
             allow_exceptions=False, label=False, id_column=None, label_file=None, label_column=None, 
             by_id=None, inference=False, protacs_embedding=False, ligand_docked=False, xyz=False):
    feature = {"node": None, "edge": None}
    seq_node_feat, seq_len = get_seq_feature(pdb_file)
    values = get_struc_feat(pdb_file=pdb_file, true_pdb=true_pdb, seq_len=seq_len, allow_exceptions=allow_exceptions,\
                                          protacs_embedding=protacs_embedding, ligand_docked=ligand_docked)
    new_pdb_file, struc_node_feat = values[0], values[-1]
    covalent = struc_node_feat['atom_emb']['covalent']
    if not label and not inference:
        true_node_feat = get_struc_feat(pdb_file=true_pdb, true_pdb=true_pdb, seq_len=seq_len, allow_exceptions=allow_exceptions,\
                                          protacs_embedding=protacs_embedding, ligand_docked=ligand_docked)[-1]
    elif not label and inference:
        label_ = label = -1
        true_node_feat = None
        true_pos = None

    else:
        df = pd.read_csv(label_file)
        label_ = df.loc[df[id_column] == by_id, label_column].values[0]

        if not isinstance(label_, float):
            label_ = int(label_)
        else:
            label_ = float(label_)
        print('label: ', label, '\n', 'label =', label_)
        true_node_feat = None
        true_pos = None

    seq_node_feat, seq_len = get_seq_feature(new_pdb_file, updated=True)
    res_num = seq_node_feat['onehot'].shape[0]
    com_emb = []
    atom_xyz = []
    print("total resnums is: ", res_num)
    print("embedded size - resnume: ", len(struc_node_feat['atom_emb']['embedding']))
    for i in range(res_num):
        seq_emb = seq_node_feat['onehot'][i]
        atom_emb = struc_node_feat['atom_emb']['embedding'][i]
        seq_emb = np.tile(seq_emb,(atom_emb.shape[0],1))
        com_emb.append(np.concatenate((seq_emb,atom_emb[:,0:-3]),axis=1))
        atom_xyz.append(atom_emb[:,-3:])
    com_emb = np.vstack(com_emb)
    atom_xyz = np.vstack(atom_xyz)
    atom_pos = struc_node_feat['atom_emb']['atom_pos']
    atom_pos = np.vstack(atom_pos)
    if not label:
        true_pos = true_node_feat['atom_emb']['atom_pos']
        true_pos = np.vstack(true_pos)
    else:
        true_pos = label_
    atom_lst = struc_node_feat['atom_emb']['atom_lst']
    res_atom = struc_node_feat['atom_emb']['res_atom']
    CA_lst = struc_node_feat['atom_emb']['CA_lst']
    CA_pos = np.vstack(CA_lst)
    print("Building lframe pairs for all atoms ...")
    p, q, k, t = set_lframe(new_pdb_file, atom_xyz, atom_lst)
    pairs = np.concatenate([p,q,k,t],axis=-1)
    if not xyz:
        return [com_emb, atom_xyz, atom_pos, atom_lst, CA_pos, true_pos, new_pdb_file, pairs, covalent]
    else:
        return [com_emb, atom_xyz, atom_pos, atom_lst, CA_pos, true_pos, new_pdb_file, pairs, covalent, res_atom]


def main():
    parser = argparse.ArgumentParser(description='Process some pdb files.')
    parser.add_argument('pdb_file', type=str, help='Input pdb file path')
    parser.add_argument('--true_pdb', type=str, help='True pdb file path')
    parser.add_argument('--ext', default='out', type=str, help='Extension for the pickle file output')
    parser.add_argument('--label', default=False, type=bool, help='Indicate label')
    parser.add_argument('--label_file', type=str, help='Path to the label file')
    parser.add_argument('--label_column', type=str, help='Name of the label column in the label file')
    parser.add_argument('--by_id', type=str, help='ID to select specific rows or data')
    parser.add_argument('--id_column', type=str, help='ID of the identifier column')
    parser.add_argument('--allow_exceptions', default=True, type=bool, help='Indicate allow_exceptions')
    parser.add_argument('--protacs_embedding', default=False, type=bool, help='Flag to indicate if there are humongous ligands to process and/or ternary complexes')
    parser.add_argument('--ligand_docked', default=False, type=bool, help='Flag to indicate if the ligands should be parsed as one entity despite resnumber')
    parser.add_argument('--inference', default=False, type=bool, help='Flag to indicate if the script is being run for inference')
    parser.add_argument('--compress', default=True, type=bool, help='Flag to indicate whether to compress embeddings')
    parser.add_argument('--interface', type=int, help='Number of pocket residues adjacent')
    parser.add_argument('--align', default=False, type=bool, help='Flag to indicate whether to APO-HOLO preprocess')
    parser.add_argument('--clean', default=False, type=bool, help='Clean up files true/false')
    parser.add_argument('--featurize', default=False, type=bool, help='Calculate esm and rdkit descriptors (includes morgan 2048 bits) true/false')
    parser.add_argument('--scaler', type=str, help='Path to features scaler')
    parser.add_argument('--ligands_path', type=str, help="Ligands path in case rdkit fails to calculate features directly from pdb input file. Expects both sdf and mol2 files to try out and expects pattern <by_id_ligand.sdf> etc for mol2")
    parser.add_argument('--big_hets', default=False, type=bool, help='Only keep the biggest ligands in the two input pdb files true/false')
    parser.add_argument('--xyz', default=False, type=bool, help='Are we predicting new XYZ eventually true/false')


    print('starting processing ...')
    args = parser.parse_args()
    overwrite_pdb_with_pymol(args.pdb_file)
    #clean_pdb(args.pdb_file)
    #args.pdb_file = f"{args.pdb_file.split('.')[0]}_clean.pdb"
    #args.pdb_file = f"{args.pdb_file.split('.')[0]}.pdb"
    if not args.interface:
        with open(args.pdb_file, 'r') as file:
            lines = file.readlines()
            if len(lines) >= 4500:
                args.interface = 400
                args.align = True
    if args.interface:
        nearest_residues_around_ligand(args.pdb_file, f"{args.pdb_file.split('.')[0]}_interface.pdb", args.interface)
        args.pdb_file = f"{args.pdb_file.split('.')[0]}_interface.pdb"
    print("Processing PDB file: ", args.pdb_file)
    if args.true_pdb: 
        clean_pdb(args.true_pdb)
        args.true_pdb = f"{args.true_pdb.split('.')[0]}_clean.pdb"
        if args.interface:
            nearest_residues_around_ligand(args.true_pdb, f"{args.true_pdb.split('.')[0]}_interface.pdb", args.interface)
            args.true_pdb = f"{args.true_pdb.split('.')[0]}_interface.pdb"
        if args.align: 
            pdb_file, true_pdb = align_and_filter(args.pdb_file, args.true_pdb,  big_hets=args.big_hets)
            dst_dir_path = os.getcwd()
            shutil.move(pdb_file, dst_dir_path)       
            shutil.move(true_pdb, dst_dir_path)    
            args.pdb_file = pdb_file.split('/')[-1]
            args.true_pdb = true_pdb.split('/')[-1] 
    
    x = get_item(pdb_file=args.pdb_file, true_pdb=args.true_pdb, label=args.label, allow_exceptions=args.allow_exceptions,\
                label_file=args.label_file, label_column=args.label_column, by_id=args.by_id, id_column=args.id_column, \
                 inference=args.inference, protacs_embedding=args.protacs_embedding, ligand_docked=args.ligand_docked, xyz=args.xyz)
    if args.ligand_docked:
        assert x[0][-1][:25][-1] == 1, "NO LIGANDS WERE EMBEDDED"

    if args.true_pdb: 
        assert len(x[0]) == len(x[1]) == len(x[2]) == len(x[7]) == len(x[8]) == len(x[5]), "Mismatch of number of atoms across embeddings"
    else:
        assert len(x[0]) == len(x[1]) == len(x[2]) == len(x[7]) == len(x[8]), "Mismatch of number of atoms across embeddings"
    assert len(x[3]) == len(x[4]), "Mismatch of number of residues across embeddings"

    #if not x[1].shape[0] < 6000:
    #    raise ValueError("Size of the first dimension should be less than 6000")

    if args.featurize:
        features = pdb_to_features(args.pdb_file, args.scaler, args.ligands_path, args.by_id)
        x.append(features)

    if isinstance(x, list) and isinstance(x[0], np.ndarray):
        # Iterate over the numpy array and find the index
        for i in range(x[0].shape[0]):
            if x[0][i][:25][-1] == 1:
                break
        else:
            # If the loop did not break, there is no such index
            i = None

        # Append the index to x
        x.append(i)


    if not args.compress: 
        dump(x, f"{os.path.dirname(os.path.abspath(args.pdb_file))}/{x[6].split('.')[0]}_{args.ext}.pkl")
    else:
        dump(x, f"{os.path.dirname(os.path.abspath(args.pdb_file))}/{x[6].split('.')[0]}_{args.ext}.pkl", compress=3)


    file_patterns = ["*clean*pdb", "*interface.pdb", "*lframe.pdb", "*unique_chains.pdb", "*ligchain.pdb"]

    if args.clean:
        for file_pattern in file_patterns:
            for file in glob.glob(file_pattern):
                if "original.pdb" in file:  # Skip deletion if it's an original.pdb file
                    continue

                try:
                    os.remove(file)
                except OSError as e:
                    print(f"Error: {file} : {e.strerror}")


if __name__ == '__main__':
    warnings.filterwarnings("ignore")
    main()
