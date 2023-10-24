"""
------------------------------------------------------------------------------------------
Copyright (C) 2023 Hazem Mslati. All Rights Reserved.

This script is provided "AS IS" without any warranty.
Unauthorized copying, modification, or distribution is prohibited.
------------------------------------------------------------------------------------------
"""

import os
from Bio.PDB import PDBIO, PDBParser
from Bio import pairwise2
from Bio.Align import PairwiseAligner
import pymol2
from Bio.PDB.Chain import Chain
from Bio.Seq import Seq
from Bio import AlignIO



# Initialize PDB Parser
parser = PDBParser(QUIET=True)

# Directory paths
ternary_dir = "ternaries"
e3_dir = "e3"

# Get all files from directories
ternary_files = os.listdir(ternary_dir)
e3_files = os.listdir(e3_dir)


def load_and_align(e3_path, ternary_path):
    # Initialize PDB Parser
    parser = PDBParser(QUIET=True)

    # Parse structures
    e3_structure = parser.get_structure('e3', e3_path)
    ternary_structure = parser.get_structure('ternary', ternary_path)

    # Get chain IDs from Biopython structures
    e3_chain = list(e3_structure[0].get_chains())[0].id
    ternary_chain = list(ternary_structure[0].get_chains())[0].id

    with pymol2.PyMOL() as pymol:
        pymol.cmd.load(e3_path, "e3")
        pymol.cmd.load(ternary_path, "ternary")

        # Align structures using chain IDs extracted by Biopython
        pymol.cmd.align(f'ternary and chain {ternary_chain}', f'e3 and chain {e3_chain}', object="aln")

        # Create selection of matched residues, excluding HETATM and considering only CA proximity
        pymol.cmd.select('ca_match', 'ternary and polymer and name CA and (br. all within 2 of e3 and name CA)')

        # Expand selection to the entire residues
        pymol.cmd.select('matched', 'byres ca_match')

        # Change chain of matched residues
        pymol.cmd.alter('matched', 'chain="B"')

        # Save the modified structure
        ternary_file_name = os.path.basename(ternary_path)
        pymol.cmd.save(os.path.join('ternaries', f"new_{ternary_file_name}"), 'ternary')

for e3_file in e3_files:
    e3_file_name, e3_file_ext = os.path.splitext(e3_file)
    if e3_file_ext.lower() != '.pdb':  # Skip non-PDB files
        continue

    var_part = "_".join(e3_file_name.split(".")[0].split("_")[2:4]).lower() # extract the identifying part

    for ternary_file in ternary_files:
        ternary_file_name, ternary_file_ext = os.path.splitext(ternary_file)
        if ternary_file_ext.lower() != '.pdb':  # Skip non-PDB files
            continue

        var_part_ternary = "_".join(ternary_file_name.split("_")[5:7]).lower()
        if var_part_ternary == var_part:  # If the identifying parts match
            load_and_align(os.path.join(e3_dir, e3_file), os.path.join(ternary_dir, ternary_file))




