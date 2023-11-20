"""
------------------------------------------------------------------------------------------
Copyright (C) 2023 Hazem Mslati. All Rights Reserved.

This script is provided "AS IS" without any warranty.
Unauthorized copying, modification, or distribution is prohibited.
------------------------------------------------------------------------------------------
"""

import argparse
import Bio.PDB
import numpy as np
import sys
import os
import shutil
from itertools import product
import pandas as pd
import math as m
import glob
import statistics
import collections
from scipy.spatial import distance
from biopandas.pdb import PandasPdb
from itertools import groupby, count
from ast import literal_eval

parser = Bio.PDB.PDBParser(QUIET=True)

pdb_docked= str(sys.argv[1])

ppdb = PandasPdb().read_pdb(pdb_docked)
lig_x_coord,lig_y_coord,lig_z_coord= statistics.mean(list(ppdb.df['HETATM'].x_coord)), statistics.mean(list(ppdb.df['HETATM'].y_coord)), statistics.mean(list(ppdb.df['HETATM'].z_coord))
ligand_coordinates_avg=[[lig_x_coord,lig_y_coord,lig_z_coord]]

protein_coord=[]
for i in range(len(ppdb.df['ATOM'])):
    protein_coord.append([ppdb.df['ATOM'].x_coord[i], ppdb.df['ATOM'].y_coord[i], ppdb.df['ATOM'].z_coord[i]])

def get_minimum_residue_distance (avg_ligand_coor, rec_coors):
    dis_df = pd.DataFrame(columns=["A", "B", "distance"])

    for pair in product(avg_ligand_coor, rec_coors):
        x, y = pair[0], pair[1]

        dist = distance.euclidean(x, y)
        dis_df = dis_df.append(
            {'A': x, 'B': y, 'distance': dist}, ignore_index=True
        )
    the_minimum_index=dis_df.index[dis_df.distance == dis_df.distance.min()]
    return(the_minimum_index.tolist()[0])

most_minimum_dist_residue_index=get_minimum_residue_distance (ligand_coordinates_avg, protein_coord)
def parse_pdb_structure (pdb):
    return(parser.get_structure(str(pdb.rsplit( ".", 1 )[ 0 ]) , pdb))

def what_chain_is_poi (structure):
    chains=[]
    for model in structure:
        for chain in model:
            chains.append(chain)
    chain_id=int()
    for i in range(len(chains)):
        residues_list=[]
        for idx, residue in enumerate(chains[i]):
            residues_list.append(residue)
        if len(residues_list)>=30:
            chain_id=i
            break
        else:
            print(f'poi not chain {i+2} ??')
            continue
    return(chain_id)

chain_id=what_chain_is_poi(parse_pdb_structure(pdb_docked))

target_residue_name=ppdb.df['ATOM'].iloc[most_minimum_dist_residue_index]['residue_name']
target_residue_number=ppdb.df['ATOM'].iloc[most_minimum_dist_residue_index]['residue_number']
target_residue=f'{target_residue_name}{target_residue_number}'

structure=parse_pdb_structure(pdb_docked)
chains=[]
idx_model=int()
for model in structure:
    for chain in model:
        chains.append(chain)

    for idx, residue in enumerate(chains[0]):
        if f'resseq={target_residue_number} ' in str(residue):
            idx_model=idx
            
residues_ref = [r for r in structure.get_residues()]
target_atom = residues_ref[idx_model]['CA']
atoms = Bio.PDB.Selection.unfold_entities(structure, 'A')
ns = Bio.PDB.NeighborSearch(atoms)
close_residues = ns.search(target_atom.coord, 10 ,'R')
# close_atoms=[coor.coord for coor in close_atoms]
res_all=[close_residues[i].get_full_id()[3][1] for i in range(len(close_residues))]

res_all_s=sorted(res_all)

def as_range(g):
    l = list(g)
    return l[0], l[-1]

res_all_p=[as_range(g) for _, g in groupby(res_all_s, key=lambda n, c=count(): n-next(c))]

res_no_dup=[sorted(list(set(p))) for p in res_all_p]

res_all_str=[]
count=0
for ele in res_no_dup:
    str1 = ""
    count=len(ele)
    for elem in ele:
        str1 += str(elem)
        while count > 1:
            str1 +='-'
            count-=1
            break
        else:
            continue
    res_all_str.append(str1)


textfile = open(str(sys.argv[2]), "w")
for element in res_all_str:
    textfile.write(element + ",")
textfile.close()

