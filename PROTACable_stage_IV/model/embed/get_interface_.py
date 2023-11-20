"""
------------------------------------------------------------------------------------------
Copyright (C) 2023 Hazem Mslati. All Rights Reserved.

This script is provided "AS IS" without any warranty.
Unauthorized copying, modification, or distribution is prohibited.
------------------------------------------------------------------------------------------
"""


from Bio.PDB import *
import argparse
from operator import itemgetter
import numpy as np
from collections import defaultdict

def calc_residue_dist(residue_one, residue_two):
    diff_vector = residue_one["CA"].coord - residue_two["CA"].coord
    return np.sqrt(np.sum(diff_vector * diff_vector))

def is_contact(residue_one, residue_two, threshold=8.0):
    return calc_residue_dist(residue_one, residue_two) < threshold

def clean_pdb(filename):
    residue_atoms = []
    last_residue_num = -1
    last_resname = ""
    lines_to_keep = []
    current_residue_index = 1
    current_atom_index = 1  # ATOM serial numbers usually start from 1
    het_residue_index = 1   # separate counter for HETATM residue numbers
    chain_id = None
    with open(filename) as f:
        for line in f:
            record_name = line[:6].strip()
            if record_name in {"ATOM", "HETATM"}:
                if record_name == "ATOM":
                    # Get chain_id from first ATOM record
                    if chain_id is None:
                        chain_id = line[21]
                    if current_residue_num != last_residue_num and last_residue_num != -1:
                        if sum(1 for atom in residue_atoms if not atom[13:14].startswith('H')) >= 3:
                            for atom_line in residue_atoms:
                                atom_line = atom_line[:6] + f"{current_atom_index:>5}" + atom_line[11:22] + f"{current_residue_index:>4}" + atom_line[26:]
                                lines_to_keep.append(atom_line)
                                current_atom_index += 1
                            current_residue_index += 1
                        residue_atoms = []
                    last_residue_num = current_residue_num
                    residue_atoms.append(line)
                else:  # record_name == "HETATM"
                    # Before we handle the HETATM, we should append the atoms from the last residue, if any.
                    if residue_atoms:
                        if sum(1 for atom in residue_atoms if not atom[13:14].startswith('H')) >= 3:
                            for atom_line in residue_atoms:
                                atom_line = atom_line[:6] + f"{current_atom_index:>5}" + atom_line[11:22] + f"{current_residue_index:>4}" + atom_line[26:]
                                lines_to_keep.append(atom_line)
                                current_atom_index += 1
                            current_residue_index += 1
                        residue_atoms = []
                    if current_resname != last_resname:
                        het_residue_index += 1
                    last_resname = current_resname
                    line = line[:6] + f"{current_atom_index:>5}" + line[11:22] + f"{current_residue_index + het_residue_index + 1:>4}" + line[26:]
                    lines_to_keep.append(line)
                    current_atom_index += 1
            elif record_name == "TER":
                # Before we handle the TER line, we should append the atoms from the last residue, if any.
                if residue_atoms:
                    if sum(1 for atom in residue_atoms if not atom[13:14].startswith('H')) >= 3:
                        for atom_line in residue_atoms:
                            atom_line = atom_line[:6] + f"{current_atom_index:>5}" + atom_line[11:22] + f"{current_residue_index:>4}" + atom_line[26:]
                            lines_to_keep.append(atom_line)
                            current_atom_index += 1
                        current_residue_index += 1
                    residue_atoms = []
                line = line[:6] + f"{current_atom_index:>5}" + line[11:]
                lines_to_keep.append(line)
                current_atom_index += 1
    lines_to_keep_updated = []
    for line in lines_to_keep:
        if line.startswith("HETATM") and line[21] == ' ':
            # The chain ID is inserted after the residue name and before the residue number
            newline = line[:21] + chain_id + line[22:]
            lines_to_keep_updated.append(newline)
        else:
            lines_to_keep_updated.append(line)

    with open(f"{filename.split('.')[0]}_clean.pdb", "w") as f:
        for line in lines_to_keep_updated:
            f.write(line)
    lines_to_keep.append("END\n")
    with open(f"{filename.split('.')[0]}_clean.pdb", "w") as f:
        for line in lines_to_keep:
            f.write(line)


def get_coords_from_line(line):
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    return x, y, z

def calc_centroid(hetatms):
    coords = np.array([get_coords_from_line(atom) for atom in hetatms])
    length = coords.shape[0]
    centroid = np.sum(coords, axis=0) / length
    return centroid

def get_hetatms(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    hetatms = defaultdict(list)
    for line in lines:
        if line.startswith('HETATM'):
            resname = line[17:20].strip()
            hetatms[resname].append(line)
    largest_hetatm_resname = max(hetatms, key=lambda k: len(hetatms[k]))
    return hetatms[largest_hetatm_resname]

def get_hetatms_(file_path):
    hetatms = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('HETATM'):
                hetatms.append(line)
    return hetatms


def calc_residue_dist(residue, ligand_centroid):
    min_dist = 1e6  # large initial minimum distance
    for atom_res in residue:
        if atom_res.name == "CA":  # only compare based on alpha-carbon
            diff_vector = atom_res.coord - ligand_centroid
            dist = np.sqrt(np.sum(diff_vector * diff_vector))
            if dist < min_dist:
                min_dist = dist
    return min_dist

def distance(coord1, coord2):
    return np.linalg.norm(np.array(coord1) - np.array(coord2))

def filter_far_hetatms(hetatms):
    # Extract coordinates of all HETATM lines
    coordinates = [get_coords_from_line(line) for line in hetatms]
    
    to_keep = []

    for i, coord in enumerate(coordinates):
        keep_current = False
        for j, other_coord in enumerate(coordinates):
            if i != j and distance(coord, other_coord) <= 7.5:
                keep_current = True
                break
        if keep_current:
            to_keep.append(hetatms[i])

    return to_keep


def nearest_residues_around_ligand(input_pdb, output_pdb, num_res):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('pdb', input_pdb)

    hetatms = get_hetatms(input_pdb)
    hetatms_ = get_hetatms_(input_pdb)
    filtered_hetatms = filter_far_hetatms(hetatms_)

    ligand_centroid = calc_centroid(hetatms)

    dist_list = []
    for model in structure:
        for chain in model:
            for residue in chain:
                dist_list.append((residue, calc_residue_dist(residue, ligand_centroid)))

    dist_list.sort(key=itemgetter(1))

    nearest_residues = [item[0] for item in dist_list[:num_res]]

    io = PDBIO()
    io.set_structure(structure)

    class ResidueSelect(Select):
        def accept_residue(self, residue):
            # Check if the residue is a ligand
            if residue.id[0].startswith('H_'):
                return 0
            # Check if the residue is in the nearest_residues list
            return residue in nearest_residues

    io.save(output_pdb, ResidueSelect())

    with open(output_pdb, 'r') as f:
        lines = f.readlines()
    with open(output_pdb, 'w') as f:
        f.writelines(lines[:-1])

    with open(output_pdb, 'a') as f:
        for line in filtered_hetatms:
            f.write(line)
        f.write('END\n')

if __name__ == "__main__":
    # Argument parsing
    parser = argparse.ArgumentParser(description='Find the nearest residues around a ligand.')
    parser.add_argument('input_pdb', type=str, help='Input PDB file path.')
    parser.add_argument('output_pdb', type=str, help='Output PDB file path.')
    parser.add_argument('num_res', type=int, help='Number of nearest residues to keep.')

    args = parser.parse_args()
    nearest_residues_around_ligand(args.input_pdb, args.output_pdb, args.num_res)

