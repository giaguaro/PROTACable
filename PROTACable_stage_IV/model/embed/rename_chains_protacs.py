"""
------------------------------------------------------------------------------------------
Copyright (C) 2023 Hazem Mslati. All Rights Reserved.

This script is provided "AS IS" without any warranty.
Unauthorized copying, modification, or distribution is prohibited.
------------------------------------------------------------------------------------------
"""

from multiprocessing import Pool
import os
from Bio.PDB import PDBIO, PDBParser
from Bio import pairwise2
from Bio.Align import PairwiseAligner
import pymol2
from Bio.PDB.Chain import Chain
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.SeqUtils import seq1
from Bio.PDB import PPBuilder
from Bio import Align
from Bio.SubsMat import MatrixInfo as matlist
from Bio import SeqIO
from Bio.PDB import *
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB import PDBIO, PDBParser, Model, Chain, Residue

# Initialize the aligner
aligner = Align.PairwiseAligner()

# Initialize the PPBuilder class
ppb = PPBuilder()

# Initialize PDB Parser
parser = PDBParser(QUIET=True)

# Directory paths
ternary_dir = "ternaries"
e3_dir = "e3"
poi_dir = "poi"

# Get all files from directories
ternary_files = os.listdir(ternary_dir)
e3_files = os.listdir(e3_dir)
poi_files = os.listdir(poi_dir)

def get_sequence(structure):
    for model in structure:
        for chain in model:
            seq = [residue.get_resname() for residue in chain]
            return seq  # return sequence from the first chain


# List to store matched files
matched_files = []

# Iterate over e3 files
for e3_file in e3_files:
    e3_file_name, e3_file_ext = os.path.splitext(e3_file)
    if e3_file_ext.lower() != '.pdb':  # Skip non-PDB files
        continue

    var_part = "_".join(e3_file_name.split(".")[0].split("_")[2:4]).lower() # extract the identifying part

    for poi_file in poi_files:
        poi_file_name, poi_file_ext = os.path.splitext(poi_file)
        if poi_file_ext.lower() != '.pdb':  # Skip non-PDB files
            continue

        var_part_poi = "_".join(poi_file_name.split(".")[0].split("_")[1:3]).lower()

        for ternary_file in ternary_files:
     #       ternary_file_name, ternary_file_ext = os.path.splitext(ternary_file)
      #      if ternary_file_ext.lower() != '.pdb':  # Skip non-PDB files
       #         continue

            var_part_ternary = ternary_file #"_".join(ternary_file_name.split("_")[5:7]).lower()

            if var_part_ternary == var_part == var_part_poi:  # If the identifying parts match
                matched_files.append([poi_file, e3_file, ternary_file])

                
def process_matched_file(matched_file):
    
    directory_path = f"{os.path.join(ternary_dir,matched_file[2])}/processed_{matched_file[2]}/"
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)

    current_ternary = os.path.join(ternary_dir,matched_file[2])
    ternary_subfiles = [file for file in os.listdir(current_ternary) if file.endswith(".pdb")]

    for pose in ternary_subfiles:

        parser = PDBParser()

        poi_structure = parser.get_structure("POI",  os.path.join(poi_dir, matched_file[0]))
        e3_structure = parser.get_structure("E3", os.path.join(e3_dir, matched_file[1]))
        ternary_structure = parser.get_structure("TERNARY",  os.path.join(current_ternary, pose))

        sequence1 = seq1("".join(get_sequence(ternary_structure)))

        sequence2 = seq1("".join(get_sequence(poi_structure)))

        sequence3 = seq1("".join(get_sequence(e3_structure)))

        alignments_poi = pairwise2.align.globalxx(sequence1, sequence2)
        alignment_poi = alignments_poi[0]

        alignments_e3 = pairwise2.align.globalxx(sequence1, sequence3)
        alignment_e3 = alignments_e3[0]

        # Find index of the first aligned pair
        all_residues = list(ternary_structure[0].get_residues())
        res_index_map = {residue.id[1]: i for i, residue in enumerate(all_residues)}

        first_aligned_pair_index = next((i for i, pair in enumerate(zip(alignment_e3[0], alignment_e3[1])) if (pair[0] != '-' and pair[1] != '-')), None)
        print(first_aligned_pair_index)
 #       first_aligned_residue = all_residues[res_index_map[first_aligned_pair_index]]
#        print(first_aligned_pair_index)

        # Check if we found an index
        if first_aligned_pair_index is not None:
            aligned_structure = Structure.Structure('ALIGNED')

            # Create a new model in the structure
            model = Model.Model(0)
            aligned_structure.add(model)

            # Create a new chain in the model
            chainA = Chain.Chain('A')
            chainB = Chain.Chain('B')
            model.add(chainA)
            model.add(chainB)

            # Add each part's residues to the respective Model
            for i, res in enumerate(all_residues):
                if i < first_aligned_pair_index:
                    try:
                        res = ternary_structure[0]['A'][i]
                        new_res = Residue.Residue(res.get_id(), res.get_resname(), res.get_segid())
                        for atom in res:
                            new_res.add(atom.copy())
                        chainA.add(new_res)
                    except:
                        pass
                elif i >= first_aligned_pair_index:
                    try:
                        res = ternary_structure[0]['A'][i]
                        new_res = Residue.Residue(res.get_id(), res.get_resname(), res.get_segid())
                        for atom in res:
                            new_res.add(atom.copy())
                        chainB.add(new_res)
                    except:
                        pass

            io = PDBIO()
            io.set_structure(aligned_structure)
            io.save(f"{directory_path}/aligned_{os.path.splitext(pose)[0]}.pdb")




# Create a multiprocessing Pool
pool = Pool()

# Create a list of arguments for each task
pool.map(process_matched_file, matched_files)

