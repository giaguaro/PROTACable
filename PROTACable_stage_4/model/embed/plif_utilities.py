"""
------------------------------------------------------------------------------------------
Copyright (C) 2023 Hazem Mslati. All Rights Reserved.

This script is provided "AS IS" without any warranty.
Unauthorized copying, modification, or distribution is prohibited.
------------------------------------------------------------------------------------------
"""

# plif_utilities.py

from tqdm import tqdm
from rdkit.Chem import rdMolTransforms, rdMolDescriptors
from rdkit.Geometry import Point3D
from rdkit import Chem
from rdkit.Chem import AllChem
import copy
import tempfile
import Bio.PDB
from Bio.PDB import PDBIO
from Bio.PDB import PDBParser, Vector, PDBIO, Select
from os import system
import openbabel
import numpy as np
from collections import defaultdict
from biopandas.pdb import PandasPdb
import math
from itertools import groupby
from sklearn.decomposition import PCA
from joblib import Parallel, delayed




SUBSTITUTIONS = {
    '2AS':'ASP', '3AH':'HIS', '5HP':'GLU', '5OW':'LYS', 'ACL':'ARG', 'AGM':'ARG', 'AIB':'ALA', 'ALM':'ALA', 'ALO':'THR', 'ALY':'LYS', 'ARM':'ARG',
    'ASA':'ASP', 'ASB':'ASP', 'ASK':'ASP', 'ASL':'ASP', 'ASQ':'ASP', 'AYA':'ALA', 'BCS':'CYS', 'BHD':'ASP', 'BMT':'THR', 'BNN':'ALA',
    'BUC':'CYS', 'BUG':'LEU', 'C5C':'CYS', 'C6C':'CYS', 'CAS':'CYS', 'CCS':'CYS', 'CEA':'CYS', 'CGU':'GLU', 'CHG':'ALA', 'CLE':'LEU', 'CME':'CYS',
    'CSD':'ALA', 'CSO':'CYS', 'CSP':'CYS', 'CSS':'CYS', 'CSW':'CYS', 'CSX':'CYS', 'CXM':'MET', 'CY1':'CYS', 'CY3':'CYS', 'CYG':'CYS',
    'CYM':'CYS', 'CYQ':'CYS', 'DAH':'PHE', 'DAL':'ALA', 'DAR':'ARG', 'DAS':'ASP', 'DCY':'CYS', 'DGL':'GLU', 'DGN':'GLN', 'DHA':'ALA',
    'DHI':'HIS', 'DIL':'ILE', 'DIV':'VAL', 'DLE':'LEU', 'DLY':'LYS', 'DNP':'ALA', 'DPN':'PHE', 'DPR':'PRO', 'DSN':'SER', 'DSP':'ASP',
    'DTH':'THR', 'DTR':'TRP', 'DTY':'TYR', 'DVA':'VAL', 'EFC':'CYS', 'FLA':'ALA', 'FME':'MET', 'GGL':'GLU', 'GL3':'GLY', 'GLZ':'GLY',
    'GMA':'GLU', 'GSC':'GLY', 'HAC':'ALA', 'HAR':'ARG', 'HIC':'HIS', 'HIP':'HIS', 'HMR':'ARG', 'HPQ':'PHE', 'HTR':'TRP', 'HYP':'PRO',
    'IAS':'ASP', 'IIL':'ILE', 'IYR':'TYR', 'KCX':'LYS', 'LLP':'LYS', 'LLY':'LYS', 'LTR':'TRP', 'LYM':'LYS', 'LYZ':'LYS', 'MAA':'ALA', 'MEN':'ASN',
    'MHS':'HIS', 'MIS':'SER', 'MK8':'LEU', 'MLE':'LEU', 'MPQ':'GLY', 'MSA':'GLY', 'MSE':'MET', 'MVA':'VAL', 'NEM':'HIS', 'NEP':'HIS', 'NLE':'LEU',
    'NLN':'LEU', 'NLP':'LEU', 'NMC':'GLY', 'OAS':'SER', 'OCS':'CYS', 'OMT':'MET', 'PAQ':'TYR', 'PCA':'GLU', 'PEC':'CYS', 'PHI':'PHE',
    'PHL':'PHE', 'PR3':'CYS', 'PRR':'ALA', 'PTR':'TYR', 'PYX':'CYS', 'SAC':'SER', 'SAR':'GLY', 'SCH':'CYS', 'SCS':'CYS', 'SCY':'CYS',
    'SEL':'SER', 'SEP':'SER', 'SET':'SER', 'SHC':'CYS', 'SHR':'LYS', 'SMC':'CYS', 'SOC':'CYS', 'STY':'TYR', 'SVA':'SER', 'TIH':'ALA',
    'TPL':'TRP', 'TPO':'THR', 'TPQ':'ALA', 'TRG':'LYS', 'TRO':'TRP', 'TYB':'TYR', 'TYI':'TYR', 'TYQ':'TYR', 'TYS':'TYR', 'TYY':'TYR',
    "5MC": "C",  # 5-methylcytosine
    "5HM": "C",  # 5-hydroxymethylcytosine
    
    # RNA modifications
    "OMC": "C",  # O-methylcytidine
    "OMG": "G",  # O-methylguanosine
    "OMU": "U",  # O-methyluridine
    "O2M": "U",  # 2'O-methylated nucleotide; could be various bases, here I'm approximating to U
    "PSU": "U",  # Pseudouridine
    "H2U": "U",  # Dihydrouridine
    "M2G": "G",  # N2-methylguanosine
    "I": "G",    # Inosine, often derived from Adenosine but more similar in structure to Guanosine
    "M1A": "A",  # 1-methyladenosine
    "M7G": "G",  # 7-methylguanosine
    "T6A": "A",  # N6-threonylcarbamoyladenosine
    "M6A": "A",  # N6-methyladenosine
    "M5C": "C",  # 5-methylcytidine
    "M5U": "U",  # 5-methyluridine
    "S4U": "U",  # 4-thiouridine
    "C5M": "C",  # 5-methylcytidine (alternative code)
}

COMMON_IONS = ['NA', 'K', 'MG', 'CA', 'CO', 'ZN', 'MN', 'FE', 'NI', 'HG', 'CD', 'CU', 'AG']


# Atom properties mapped to values
backbone_atom_properties = {
    0: [3, 0, 0, 0, 0],  # N
    1: [3, 0, 0, 0, 0],  # CA
    2: [2, 0, 0, 0, 0],  # C
    3: [2, 0, 0, 0, 0],  # O
}

# {atomic number : atom symbol}
atom_symbols = {1: 'H', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S',
    17: 'Cl', 19: 'K', 20: 'Ca', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn', 31: 'Ga', 33: 'As',
    35: 'Br', 42: 'Mo', 44: 'Ru', 53: 'I', 74: 'W', 79: 'Au'
}
# {atom symbol + possible hybridization state : one-hot encoding}
atom_dict_protein = {"N":0, "CA":1, "C":2, "O":3, "CB":4, "OG":5, "CG":6, "CD1":7, "CD2":8, "CE1":9, "CE2":10, "CZ":11, "OD1":12,
                     "ND2":13, "CG1":14, "CG2":15, "CD":16, "CE":17, "NZ":18, "OD2":19, "OE1":20, "NE2":21, "OE2":22, "OH":23,
                     "NE":24, "NH1":25, "NH2":26, "OG1":27, "SD":28, "ND1":29, "SG":30, "NE1":31, "CE3":32, "CZ2":33, "CZ3":34, "CH2":35,
                     "OXT":36}

backbone_atoms = ['N', 'CA', 'C', 'O']

atom_dict_ligand = {'H': 0, 'C': 1, 'N': 2, 'O': 3, 'F': 4, 'Na': 5,
                    'Mg': 6, 'Al': 7, 'Si': 8, 'P': 9, 'S': 10, 'Cl': 11, 'K': 12, 'Ca': 13,
                    'Mn': 14, 'Fe': 15, 'Co': 16, 'Ni': 17, 'Cu': 18, 'Zn': 19, 'Ga': 20, 'As': 21,
                    'Br': 22, 'Mo': 23, 'Ru': 24, 'I': 25, 'W': 26, 'Au': 27, 'Ag': 28, 'Se': 29, 'Cd': 30, 'B':31}

hybridization = {'S': 0, 'SP': 1, 'SP2': 2, 'SP3': 3, 'SP3D': 4, 'SP3D2': 5}
partial_charge = {-3: 0, -2: 1, -1: 2, 0: 3, 1: 4, 2: 5, 3: 6}
ring_membership = {'False': 0, 'True': 1}
aromaticity = {'False': 0, 'True': 1}
chirality = {'0': 0, 'R': 1, 'S': 2}
degree = {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6}
total_num_Hs = {0: 0, 1: 1, 2: 2, 3: 3, 4: 4}
implicit_valence = {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5}
formal_charge = {-4: 0, -3: 1, -2: 2, -1: 3, 0: 4, 1: 5, 2: 6, 3: 7, 4: 8}
heavy_atom_count = {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6}
in_ring_size = {0: 0, 3: 1, 4: 2, 5: 3, 6: 4, 7: 5, 8: 6}
is_in_ring = {'False': 0, 'True': 1}

# XXX is for other residues (not present in the list). LIG is for ligand molecule as a residue.
RESIDUE_TYPES = [
    'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU',
    'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR',
    'XXX', 'ION', 'LIG', 'PROTEIN', 'LIGAND'
]

DNA_RNA = ['DA', 'DC', 'DG', 'DT', 'A', 'C', 'G', 'U']
HETS = ['DA', 'DC', 'DG', 'DT', 'A', 'C', 'G', 'U', 'ION', 'LIG']

property_dict_lengths = {'hybridization': 6, 'partial_charge': 7, 'ring_membership': 2,
                         'aromaticity': 2, 'chirality': 3, 'degree': 7,
                         'total_num_Hs': 5, 'implicit_valence': 6,
                         'formal_charge': 9,
                         'heavy_atom_count': 7, 'in_ring_size': 7,
                         'is_in_ring': 2}


pdb_parser = Bio.PDB.PDBParser(QUIET = True)

def get_residues_from_pdb(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure('my_structure', pdb_file)

    res_lst = []
    for model in structure:
        for chain in model:
            for residue in chain:
                # Exclude water residues
                if residue.get_resname() != 'HOH':
                    res_lst.append(residue.get_resname())

    return res_lst

# check if the molecule is just a hydrogen atom
def is_hydrogen(mol):
    return mol.GetNumAtoms() == 1 and mol.GetAtomWithIdx(0).GetAtomicNum() == 1

def euclidean_distance(point1, point2):
    return np.linalg.norm(np.array(point1) - np.array(point2))

def sum_third_elements(arrays):
    total = 0
    for array in arrays:
        if len(array) >= 3:
            total += array[2]
        else:
            total += 1
    return total

def get_atom_emb(side_chains_rdmols, ligands_rdmols, fragment_counts, pdb_file, res_range, allow_exceptions=False):
    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure('protein', pdb_file)
    model = structure[0]

    res_lst = get_residues_from_pdb(pdb_file)
    res_num = res_range[1] - res_range[0] + 1

    if fragment_counts: total_frags = sum_third_elements(fragment_counts)
    else: total_frags = 0

    print(f"Total protein residues = {res_num}; Total ligand residues upon fragmentation = {total_frags}")

    protein_atom_embs = [-1 for _ in range(res_num + total_frags)]
    protein_atom_xyz = [-1 for _ in range(res_num + total_frags)]
    ligand_atom_embs = []
    ligand_atom_xyz = []

    protein_atom_lst = []
    ligand_atom_lst = []
    protein_res_atom = defaultdict(list)
    ligand_res_atom = defaultdict(list)
    CA_lst = []

    exclude_res = []

    # here i have excluded hydrogens to simplify things..
    if allow_exceptions:
        print("Allowing edge cases in the protein / ligand ")
        total_protein_atoms = sum(len([atom for atom in rdmol.GetAtoms() if atom.GetAtomicNum() != 1])
                                  for rdmols in side_chains_rdmols.values() for rdmol in rdmols if rdmol is not None)


        flattened_molecules = [rdmol for rdmols in side_chains_rdmols.values() for rdmol in rdmols]
        # some mols are None
        none_indices = [i for i, rdmol in enumerate(flattened_molecules) if rdmol is None]
        exclude_res = none_indices
        sum_nones = len(none_indices)
        res_num -= sum_nones

        # Convert none_indices to a set for faster lookup
        none_indices_set = set(none_indices)

        # Create a new list that excludes elements at the indices in none_indices
        res_lst = [val for i, val in enumerate(res_lst) if i not in none_indices_set]
        protein_atom_embs = protein_atom_embs[sum_nones:]
        protein_atom_xyz = protein_atom_xyz[sum_nones:]
        print("Total invalid protein residues : ", sum_nones)
    else:
        total_protein_atoms = sum(len([atom for atom in rdmol.GetAtoms() if atom.GetAtomicNum() != 1])
                                  for rdmols in side_chains_rdmols.values() for rdmol in rdmols)

    protein_covalent_bond = np.zeros((total_protein_atoms, total_protein_atoms))

    prev_side_chain_rdmol = None
    prev_atom_count = 0
    atom_counter = 0
    i_counter = 0
    chain_break = False

    for chain_id, side_chain_rdmol_lst in side_chains_rdmols.items():

        for i, side_chain_rdmol in enumerate(side_chain_rdmol_lst):
            # collateral to the exception mentioned above
            if not side_chain_rdmol: chain_break = True; print(f'Chain break detected between residues {i-1} and {i} due to invalid residue {i}.');  continue

            CA_atom = side_chain_rdmol.GetAtoms()[1] if (side_chain_rdmol.GetAtoms()[1].GetSymbol() == 'C' and side_chain_rdmol.GetAtoms()[1].GetDegree() > 1) else None
            if not CA_atom:
                for atom in side_chain_rdmol.GetAtoms():
                    if atom.GetSymbol() == 'C' and atom.GetDegree() > 1: # CA atoms usually form more than one bond
                        ref_pos = np.array(side_chain_rdmol.GetConformer().GetAtomPosition(atom.GetIdx()))
                        break
            else:
                ref_pos = np.array(side_chain_rdmol.GetConformer().GetAtomPosition(CA_atom.GetIdx()))
            CA_lst.append(ref_pos.astype(np.float16))

            if i > 0 and not chain_break:  # we can only check the distance from the second residue onwards
                try: 
                    prev_CA = CA_lst[i-1]
                    curr_CA = ref_pos
                    distance = np.linalg.norm(prev_CA - curr_CA)
 
                    if distance > 4.5:  # change this to your preferred threshold
                        chain_break = True
                        print(f'Chain break detected between residues {i-1} and {i}. Distance: {distance}Å')
                    else:
                        chain_break = False
                except:
                    chain_break = True
             

            atom_count = 0
            atom_pos, onehot, props = [], [], []

            # Initialize a counter
            atom_index = 0

            # Dictionary to store the mapping from old indices to new indices
            new_indices = {}

            for atom in side_chain_rdmol.GetAtoms():
                if atom.GetAtomicNum() == 1:
                    continue


                # Add the mapping from the old index to the new index
                new_indices[atom.GetIdx()] = atom_index
                atom_pos.append(np.array(side_chain_rdmol.GetConformer().GetAtomPosition(atom.GetIdx()))) 
                try:
                    chirality_value = str(next((chirality for atom_index, chirality in Chem.FindMolChiralCenters(side_chain_rdmol) if atom_index == atom.GetIdx()), '0'))
                except Exception as e:
                    chirality_value = '0'  # Default value in case of any error
                    print(f"Error while finding chiral centers: {e}")  # You can log or print the error if needed

                atom_properties = [
                    hybridization[atom.GetHybridization().name],
                    partial_charge[atom.GetFormalCharge()],
                    ring_membership[str(atom.IsInRing())],
                    aromaticity[str(atom.GetIsAromatic())],
                    chirality[chirality_value],
                    degree[atom.GetDegree()],
                    total_num_Hs[atom.GetTotalNumHs()],
                    implicit_valence[atom.GetImplicitValence()],
                    formal_charge[atom.GetFormalCharge()], # Formal charge
                    heavy_atom_count[sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() > 1)], # Number of attached heavy atoms
                    in_ring_size[next((size for size in range(3, 9) if atom.IsInRingSize(size)), 0)], # In ring size
                    is_in_ring[str(atom.IsInRing())] # Is in ring
                ]


                _onehot = np.zeros(len(atom_dict_ligand))
                _onehot[atom_dict_ligand[atom.GetSymbol()]] = 1
                onehot.append(np.array(_onehot))
                ohe_properties=[]
                # Create new one-hot encodings for each property
                for j, prop in enumerate(atom_properties):
                    length = list(property_dict_lengths.values())[j]
                    one_hot = np.zeros(length)
                    one_hot[int(prop)] = 1
                    ohe_properties.extend(one_hot)
                props.append(np.array(ohe_properties))
                atom_in_biopython = list(model.get_atoms())[atom_count]
                protein_res_atom[str(i)].append(atom_in_biopython.get_name().strip())

                # Increment the counter
                atom_index += 1
                atom_count += 1


            protein_atom_lst.append(atom_count)
            _atom_emb = np.concatenate([onehot, props, np.array(atom_pos - ref_pos[None,:])], axis=1)
            protein_atom_embs[i_counter] = _atom_emb.astype(np.float16)
            protein_atom_xyz[i_counter] = np.array(atom_pos).astype(np.float16)
            i_counter += 1

            # bond calculation
    #         idx1_global = atom_counter # The starting index of this side chain in the global atom count
            for bond in side_chain_rdmol.GetBonds():
                if bond.GetBeginAtom().GetAtomicNum() == 1 or bond.GetEndAtom().GetAtomicNum() == 1:
                    continue

                # Get the new indices
                idx1 = new_indices[bond.GetBeginAtomIdx()]
                idx2 = new_indices[bond.GetEndAtomIdx()]

                protein_covalent_bond[atom_counter+idx1][atom_counter+idx2] = 1
                protein_covalent_bond[atom_counter+idx2][atom_counter+idx1] = 1

            # peptide bond calculation
            if i > 0 and not chain_break:  # only calculate if the chain ID hasn't changed
                prev_chain_C_atom_idx = atom_counter - prev_atom_count + 2 # The 2nd index of the previous side chain is assumed to be backbone C
                curr_chain_N_atom_idx = atom_counter # The first atom of the current side chain is assumed to be backbone N
                assert prev_side_chain_rdmol.GetAtomWithIdx(2).GetSymbol() == 'C', f"Expected C atom, but got different atom for residue {i}"
                assert side_chain_rdmol.GetAtomWithIdx(0).GetSymbol() == 'N', f"Expected N atom, but got different atom for residue {i}"
                # this is equivalent to the above dont worry...
                protein_covalent_bond[prev_chain_C_atom_idx, curr_chain_N_atom_idx] = 1
                protein_covalent_bond[curr_chain_N_atom_idx, prev_chain_C_atom_idx] = 1


            prev_side_chain_rdmol = copy.deepcopy(side_chain_rdmol)
            prev_atom_count = atom_count
            chain_break = False

            atom_counter += atom_count

    if fragment_counts:
        atom_ligand_counter = 0 # New counter for ligand atoms
        total_ligand_atoms = sum([len([atom for atom in rdmol.GetAtoms() if atom.GetSymbol() != 'H']) for rdmol in ligands_rdmols.values() if (ligands_rdmols and rdmol)])
        ligand_covalent_bond = np.zeros((total_ligand_atoms, total_ligand_atoms))
        fragments_idx = [None] * total_frags

        for fragment_group, (ligand_name, ligand_mol) in zip(fragment_counts, ligands_rdmols.items()):
    #         if ligand_name not in res_lst: continue
            # Mapping atom indices from each fragment (_ligand_mol) to the original molecule (ligand_mol)
            atom_mapping = {}
            for i, (_ligand_name, _ligand_mol) in enumerate(fragment_group[0].items()):
                if is_hydrogen(_ligand_mol): continue 
    #             CA_neighbors = [pos for pos in CA_lst if any(euclidean_distance(pos, _ligand_mol.GetConformer().GetAtomPosition(atom.GetIdx())) <= 10 for atom in _ligand_mol.GetAtoms() if atom.GetSymbol() != 'H')]
    #             ref_pos = np.mean(CA_neighbors, axis=0)
                conf = _ligand_mol.GetConformer()
                coords = [conf.GetAtomPosition(atom.GetIdx()) for atom in _ligand_mol.GetAtoms()]
                # Convert RDKit's Point3D objects to numpy arrays
                coords = np.array([[point.x, point.y, point.z] for point in coords])
                # Calculate the centroid
                ref_pos = np.mean(coords, axis=0)
                CA_lst.append(ref_pos.astype(np.float16))
                atom_count = 0 # Reset for each ligand
                atom_pos, onehot, props = [], [], [] # Reset for each ligand

                for atom in _ligand_mol.GetAtoms():
                    match_found = False  # Track if a match is found
                    if atom.GetSymbol() != 'H':  # Consider only non-hydrogen atoms
                        _pos = _ligand_mol.GetConformer().GetAtomPosition(atom.GetIdx())

                       # non_hydrogen_atoms = [atom for atom in _ligand_mol.GetAtoms() if atom.GetSymbol() != 'H']

                        for orig_atom in ligand_mol.GetAtoms():
                            if orig_atom.GetSymbol() != 'H':  # Consider only non-hydrogen atoms
                                orig_pos = ligand_mol.GetConformer().GetAtomPosition(orig_atom.GetIdx())
                                if euclidean_distance(_pos, orig_pos) < 1e-3 and atom.GetSymbol() == orig_atom.GetSymbol():  # Threshold for position comparison
                                    atom_mapping[orig_atom.GetIdx()] = atom.GetIdx()
                                    match_found = True  # Update the flag
                                    break
                        if match_found:  # If a match is found
                            atom_pos.append(np.array(ligand_mol.GetConformer().GetAtomPosition(orig_atom.GetIdx())))

                            try:
                                chirality_value = str(next((chirality for atom_index, chirality in Chem.FindMolChiralCenters(side_chain_rdmol) if atom_index == atom.GetIdx()), '0'))
                            except Exception as e:
                                chirality_value = '0'  # Default value in case of any error
                                print(f"Error while finding chiral centers: {e}")  # You can log or print the error if needed

                            atom_properties = [
                                hybridization[atom.GetHybridization().name],
                                partial_charge[atom.GetFormalCharge()],
                                ring_membership[str(atom.IsInRing())],
                                aromaticity[str(atom.GetIsAromatic())],
                                chirality[chirality_value],
                                degree[atom.GetDegree()],
                                total_num_Hs[atom.GetTotalNumHs()],
                                implicit_valence[atom.GetImplicitValence()],
                                formal_charge[atom.GetFormalCharge()], # Formal charge
                                heavy_atom_count[sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() > 1)], # Number of attached heavy atoms
                                in_ring_size[next((size for size in range(3, 9) if atom.IsInRingSize(size)), 0)], # In ring size
                                is_in_ring[str(atom.IsInRing())] # Is in ring
                            ]

                            _onehot = np.zeros(len(atom_dict_ligand))
                            _onehot[atom_dict_ligand[orig_atom.GetSymbol()]] = 1
                            onehot.append(np.array(_onehot))
                            ohe_properties=[]
                            # Create new one-hot encodings for each property
                            for j, prop in enumerate(atom_properties):
                                length = list(property_dict_lengths.values())[j]
                                one_hot = np.zeros(length)
                                one_hot[int(prop)] = 1
                                ohe_properties.extend(one_hot)
                            props.append(np.array(ohe_properties))

                            atom_count += 1
                            atom_counter += 1

                ligand_atom_lst.append(atom_count)
                _atom_emb = np.concatenate([onehot, props, np.array(atom_pos - ref_pos[None,:])], axis=1)
                protein_atom_embs[i_counter] = _atom_emb.astype(np.float16)
                protein_atom_xyz[i_counter] = np.array(atom_pos).astype(np.float16)
                i_counter += 1

                idx1_global = atom_ligand_counter  # The starting index of this ligand in the global atom count for ligands

                for bond in _ligand_mol.GetBonds():
                    atom1, atom2 = bond.GetBeginAtom(), bond.GetEndAtom()
                    if atom1.GetSymbol() != 'H' and atom2.GetSymbol() != 'H':
                        idx1, idx2 = atom1.GetIdx(), atom2.GetIdx()
                        idx1 += idx1_global # Convert to global indices
                        idx2 += idx1_global
                        ligand_covalent_bond[idx1, idx2] = 1
                        ligand_covalent_bond[idx2, idx1] = 1

                atom_ligand_counter += atom_count
                fragments_idx[i] = idx1_global

            for fragment_bonds in fragment_group[1]:
                original_bonds = fragment_bonds[-1]
                fragment1_idx = fragment_bonds[0]
                fragment2_idx = fragment_bonds[1]
                idx1, idx2 = original_bonds
                atom1 = ligand_mol.GetAtomWithIdx(idx1)
                atom2 = ligand_mol.GetAtomWithIdx(idx2)
                if atom1.GetSymbol() != 'H' and atom2.GetSymbol() != 'H':
                    idx1 = atom_mapping[idx1] + fragments_idx[fragment1_idx]  # Convert to global indices
                    idx2 = atom_mapping[idx2] + fragments_idx[fragment2_idx]
                    ligand_covalent_bond[idx1, idx2] = 1
                    ligand_covalent_bond[idx2, idx1] = 1

    else:
        ligand_covalent_bond, total_ligand_atoms =None, None
    return protein_atom_embs, protein_atom_xyz, ligand_atom_embs, ligand_atom_xyz, protein_atom_lst, \
            ligand_atom_lst, CA_lst, protein_res_atom, ligand_res_atom, protein_covalent_bond, \
            ligand_covalent_bond, (total_protein_atoms, total_ligand_atoms), exclude_res      

def pdb_to_fragments(pdb_file, ligand_docked=False):
    with open(pdb_file, 'r') as file:
        lines = file.readlines()

    backbone_pdbs = []
    side_chains_pdbs = {}
    ligands_pdbs = {}

    backbone_pdb = ''
    side_chain_pdb = ''
    ligand_pdb = ''
    ligand_name = ''
    current_residue_number = None
    current_chain = ''

    for line in lines:
        if line.startswith('ATOM'):
            atom_type = line[12:16].strip()
            residue_type = line[17:20].strip()
            residue_number = int(line[22:26])
            chain = line[21]

            if current_residue_number != residue_number or current_chain != chain:
                if current_residue_number is not None and current_chain != '':
                    backbone_pdbs.append(backbone_pdb)
                    if current_chain not in side_chains_pdbs:
                        print('Chain: ', current_chain, 'detected')
                        side_chains_pdbs[current_chain] = []
                    side_chains_pdbs[current_chain].append(side_chain_pdb)
                side_chain_pdb = ''
                backbone_pdb = ''
                current_residue_number = residue_number
                current_chain = chain
            side_chain_pdb += line

        elif line.startswith('HETATM'):
            new_ligand_name = line[17:20].strip()
            ligand_residue_number = int(line[22:26])  # get ligand residue number
            
            if ligand_name != new_ligand_name and ligand_residue_number != current_residue_number and not ligand_docked:
                if ligand_name != '':
                    key = f"{ligand_name}_{current_residue_number}"  # create unique key for each ligand
                    if key not in ligands_pdbs:
                        ligands_pdbs[key] = ''
                    ligands_pdbs[key] += ligand_pdb

                ligand_pdb = ''
                ligand_name = new_ligand_name
                current_residue_number = ligand_residue_number
                
            elif ligand_name != new_ligand_name:
                if ligand_name != '':
                    key = f"{ligand_name}_{current_residue_number}"  # create unique key for each ligand
                    if key not in ligands_pdbs:
                        ligands_pdbs[key] = ''
                    ligands_pdbs[key] += ligand_pdb

                ligand_pdb = ''
                ligand_name = new_ligand_name
                current_residue_number = ligand_residue_number

            ligand_pdb += line

    if current_chain not in side_chains_pdbs:
        side_chains_pdbs[current_chain] = []
    side_chains_pdbs[current_chain].append(side_chain_pdb)
    key = f"{ligand_name}_{current_residue_number}"
    if key not in ligands_pdbs:
        ligands_pdbs[key] = ''
    ligands_pdbs[key] += ligand_pdb

    return side_chains_pdbs, ligands_pdbs

def merge_pdb_and_hetatm(pdb_file, hetatm_file, output_file, exclude_residue_indices=[]):
    # Read the input PDB file
    ppdb = PandasPdb().read_pdb(pdb_file)

    if exclude_residue_indices:
        # Assign a unique index to each residue across all chains
        ppdb.df['ATOM']['residue_index'] = ppdb.df['ATOM'].groupby(['chain_id', 'residue_number']).ngroup()

        # Remove the residues at the specified indices
        ppdb.df['ATOM'] = ppdb.df['ATOM'].loc[~ppdb.df['ATOM']['residue_index'].isin(exclude_residue_indices)]
        ppdb.df['ATOM'].drop(columns='residue_index', inplace=True)  # remove the temporary 'residue_index' column

    # Write the ATOM records to the output file
    ppdb.to_pdb(path=output_file,
                records=['ATOM'],
                gz=False,
                append_newline=True)

    # Get the serial number, resname, and chain identifier of the last ATOM record
    last_atom_record = ppdb.df['ATOM'].iloc[-1]
    serial_number = last_atom_record['atom_number']
    resname = last_atom_record['residue_name']
    last_chain_id = last_atom_record['chain_id']

    # Create the TER record line
    ter_line = f"TER    {serial_number+1}      {resname}\n"

    # Append the TER record to the output file
    with open(output_file, 'a') as f:
        f.write(ter_line)

    if hetatm_file is not None:
        # Renumber HETATM records and append them to the output file
        c = serial_number + 2
        r = int(last_atom_record['residue_number']) + 1
        prev_resname = None
        with open(hetatm_file, 'r') as infile, open(output_file, 'a') as outfile:
            for line in infile:
                if line.startswith("HETATM"):
                    element = line[75:78].strip()  # Extract the atom's element symbol
                    if element == 'H':  # If it's a hydrogen atom, skip this line
                        continue
                    current_resname = line[17:20]  # Extract the resname from the line
                    if current_resname != prev_resname:  # If resname has changed, increase the residue number
                        r += 1
                    line = line[:6] + f"{c:5d}" + line[11:21] + last_chain_id + f"{r:4d}" + line[26:]
                    c += 1
                    prev_resname = current_resname  # Update the previous resname for the next iteration
                outfile.write(line)

    # Append 'END' to the end of the PDB file
    with open(output_file, 'a') as f:
        f.write('END')

    return output_file

def align_atom_type(line, atom_type):
    """Ensures that the atom type is properly aligned at column 77 in a PDB record."""
    return line[:76] + f" {atom_type:>1}" + line[78:]


def clean_pdb(filename, overwrite=False):
    residue_atoms = []
    last_residue_num = -1
    last_resname = ""
    lines_to_keep = []
    current_residue_index = 1
    current_atom_index = 1
    het_residue_index = 1
    last_record_was_hetatm = False
    seen_coordinates = set()
    seen_atom_types = defaultdict(int)
    seen_atom_names = defaultdict(set)
    seen_ions = defaultdict(int)
    last_chain_id = ""

    ion_lines = []  # Store lines related to COMMON_IONS
    other_hetatm_lines = []  # Store other HETATM lines

    with open(filename) as f:
        for line in f:
            record_name = line[:6].strip()
            if record_name in {"ATOM", "HETATM"}:
                atom_name = line[12:16].strip()
                current_residue_num = int(line[22:26].strip())
                if current_residue_num < 0:
                    continue
                current_resname = line[17:20].strip()
                if current_resname == "ACE":
                    continue
                coordinates = (float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()))
                if coordinates in seen_coordinates:
                    continue
                seen_coordinates.add(coordinates)
                if record_name == "ATOM":
                    atom_type = line[76:79].strip()  # Extract atom type
                    line = align_atom_type(line, atom_type.upper())  # Correctly align atom type
                    last_chain_id = line[21]
                    last_record_was_hetatm = False
                    if current_residue_num != last_residue_num and last_residue_num != -1:
                        if sum(1 for atom in residue_atoms if not atom[13:14].startswith('H')) >= 3:
                            for atom_line in residue_atoms:
                                atom_line = atom_line[:6] + f"{current_atom_index:>5}" + atom_line[11:22] + f"{current_residue_index:>4}" + atom_line[26:]
                                lines_to_keep.append(atom_line.rstrip() + '\n')
                                current_atom_index += 1
                            current_residue_index += 1
                        residue_atoms = []
                    last_residue_num = current_residue_num
                    residue_atoms.append(line)
                    het_residue_index = current_residue_index + 1
                else:  # record_name == "HETATM"
                    atom_type = line[76:79].strip()  # Extract atom type
                    line = align_atom_type(line, atom_type.upper())  # Correctly align atom type
                    last_record_was_hetatm = True
                    if residue_atoms:
                        if sum(1 for atom in residue_atoms if not atom[13:14].startswith('H')) >= 3:
                            for atom_line in residue_atoms:
                                atom_line = atom_line[:6] + f"{current_atom_index:>5}" + atom_line[11:22] + f"{het_residue_index:>4}" + atom_line[26:]
                                lines_to_keep.append(atom_line.rstrip() + '\n')
                                current_atom_index += 1
                            het_residue_index += 1
                        residue_atoms = []
                    last_resname = current_resname
                    last_residue_num = current_residue_num
                    atom_type = line[76:78].strip()  # Extract atom type from last column

                    # Modified renaming logic
                    seen_atom_types[atom_type + current_resname] += 1
                    new_atom_name = atom_type.upper() + str(seen_atom_types[atom_type + current_resname])

                    line = line[:12] + f" {new_atom_name:<3}" + line[16:21] + last_chain_id + f"{het_residue_index:>4}" + line[26:]
                    if current_resname in COMMON_IONS:
                        ion_lines.append(line.rstrip() + '\n')
                    else:
                        other_hetatm_lines.append(line.rstrip() + '\n')
                    current_atom_index += 1
            elif record_name == "TER":
                if not last_record_was_hetatm:
                    if residue_atoms:
                        if sum(1 for atom in residue_atoms if not atom[13:14].startswith('H')) >= 3:
                            for atom_line in residue_atoms:
                                atom_line = atom_line[:6] + f"{current_atom_index:>5}" + atom_line[11:22] + f"{current_residue_index:>4}" + atom_line[26:]
                                lines_to_keep.append(atom_line.rstrip() + '\n')
                                current_atom_index += 1
                            current_residue_index += 1
                        residue_atoms = []
                    line = line[:6] + f"{current_atom_index:>5}" + line[11:]
                    lines_to_keep.append(line.rstrip() + '\n')
                    current_atom_index += 1
                seen_atom_types = defaultdict(int)  # Reset seen_atom_types dictionary for the next residue

    # When all lines are processed, first append ion_lines and then other_hetatm_lines to lines_to_keep
    lines_to_keep.extend(ion_lines)
    lines_to_keep.extend(other_hetatm_lines)
    lines_to_keep.append("END\n")

    output_filename = filename if overwrite else f"{filename.split('.')[0]}_clean.pdb"
    with open(output_filename, "w") as f:
        for line in lines_to_keep:
            f.write(line)

def process_residue_atoms(residue_atoms, current_atom_index, current_residue_index):
    if residue_atoms:
        for atom_line in residue_atoms:
            atom_line = atom_line[:6] + f"{current_atom_index:>5}" + atom_line[11:22] + f"{current_residue_index:>4}" + atom_line[26:]
            lines_to_keep.append(atom_line.rstrip() + '\n')
            current_atom_index += 1
        current_residue_index += 1
    residue_atoms.clear()
    return residue_atoms, current_atom_index, current_residue_index

def calculate_cb(N, CA, C):
    # Define the bond length, bond angles in radians
    CC_bond_length = 1.52
    NCA_bond_angle = np.radians(110.1)
    HNCA_dihedral_angle = np.radians(123.0)  # For most common trans isomers

    # Create the vectors along the bonds
    CA_N = np.array(N) - np.array(CA)
    CA_C = np.array(C) - np.array(CA)

    # Calculate the coordinates
    BC = (-1 * np.cos(NCA_bond_angle)) * CA_N + np.cos(HNCA_dihedral_angle) * CA_C

    # Normalize the vector and multiply by bond length
    norm = np.linalg.norm(BC)
    BC = (CC_bond_length / norm) * BC

    CB = BC + np.array(CA)

    return CB

class ProteinSelect(Select):
    def accept_residue(self, residue):
        return 1 if residue.id[0] == " " else 0  # return True for proteins (residue.id[0] == " ")


def get_principal_axes(atoms):
    """
    Get the principal axes of `atoms`.
    """
    if len(atoms) == 1:
        # Special case: If we only have one atom, we can't really define axes
        # based on the data. We could return standard axes, or handle this as
        # an error, depending on your needs.
        return np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])
    
    coords = np.array([atom.coord for atom in atoms])
    try:
        pca = PCA(n_components=3)
        pca.fit(coords)
    except:
        pca = PCA(n_components=2)
        pca.fit(coords)

    x_axis = pca.components_[0]
    y_axis = pca.components_[1]
    z_axis = np.cross(x_axis, y_axis)  # Compute the cross product to ensure orthogonality

    return x_axis, y_axis, z_axis


def calculate_c(N, CA):
    # Distances are in Angstroms, angles are in degrees
    avg_Ca_C_bond_length = 1.52
    avg_Ca_C_bond_angle = 109.5
    avg_N_Ca_C_bond_angle = 110.6

    # Conversion from degrees to radians for trigonometric functions
    avg_Ca_C_bond_angle = np.deg2rad(avg_Ca_C_bond_angle)
    avg_N_Ca_C_bond_angle = np.deg2rad(avg_N_Ca_C_bond_angle)

    # Vector from N to Ca
    N_Ca_vector = CA - N
    N_Ca_vector /= np.linalg.norm(N_Ca_vector)

    # Assuming an idealized peptide bond, the C atom will be in the same plane as the N and Ca atoms,
    # and it will be positioned such that the N-Ca-C bond angle is 110.6 degrees.
    # We calculate its position relative to the Ca atom and then translate it to absolute coordinates.
    C_relative = np.array([
        -np.cos(avg_N_Ca_C_bond_angle),
        np.sin(avg_N_Ca_C_bond_angle) * np.sin(avg_Ca_C_bond_angle),
        np.sin(avg_N_Ca_C_bond_angle) * np.cos(avg_Ca_C_bond_angle)
    ])
    C_relative *= avg_Ca_C_bond_length

    # Rotate C_relative into the same orientation as N_Ca_vector using Rodrigues' rotation formula
    k = np.array([0, 0, 1])  # rotation axis
    theta = np.arccos(np.dot(N_Ca_vector, k))  # angle between N_Ca_vector and k

    # Rodrigues' rotation formula
    C_relative_rotated = C_relative * np.cos(theta) + np.cross(k, C_relative) * np.sin(theta) + k * np.dot(k, C_relative) * (1 - np.cos(theta))

    # Translate the relative coordinates to absolute coordinates
    C = CA + C_relative_rotated

    return C

def set_lframe(structure_file, atom_xyz, atom_lst):
    '''
    Args:
        structure_file (string): the path of pdb structure file.
        res_range [int, int]: the start and end residue index, e.g. [1, 100].
    '''
    # load pdb file
    structure = pdb_parser.get_structure("tmp_stru", structure_file)
    all_residues = list(structure.get_residues())
    protein_residues = [res for res in all_residues if res.id[0] == ' ']
    protein_len = len(protein_residues)
    hetatm_residues = [res for res in all_residues if res.id[0] != ' ']


    # set local frame for each residue in pdb
    pdict = dict()
    pdict['N'] = []
    pdict['Ca'] = []
    pdict['C'] = []
    pdict['Cb'] = []
    pdict['x'] = []
    pdict['y'] = []
    pdict['z'] = []

    for residue in protein_residues:
        # Process protein residues
        pdict['N'].append(residue['N'].coord)
        pdict['Ca'].append(residue['CA'].coord)
        try:
            pdict['C'].append(residue['C'].coord)
        except KeyError:
            # Reconstruct C atom's position from N and CA atoms
            N = residue['N'].coord
            CA = residue['CA'].coord
            C = calculate_c(N, CA)  # Assuming a calculate_c function exists
            pdict['C'].append(C)

        # If no CB atom, calculate it
        if 'CB' in residue:
            pdict['Cb'].append(residue['CB'].coord)
        else:
            N = residue['N'].coord
            CA = residue['CA'].coord
            C = pdict['C'][-1]  # Use the most recently appended C atom
            CB = calculate_cb(N, CA, C)
            pdict['Cb'].append(CB)

    # sort hetatm_residues by resname and group them
    hetatm_residues.sort(key=lambda res: res.resname)
    hetatm_groups = {k: list(v) for k, v in groupby(hetatm_residues, key=lambda res: res.resname)}

    for hetatm_resname, hetatm_residues in hetatm_groups.items():
        hetatm_mol = [atom for res in hetatm_residues for atom in res]

        # Compute geometric centroid of the interacting atoms as the origin
        centroid = np.mean([atom.coord for atom in hetatm_mol], axis=0)
        pdict['Cb'].append(centroid)

        # Define axes based on the principal axes of the interacting atoms
        x_axis, y_axis, z_axis = get_principal_axes(hetatm_mol)
        pdict['x'].append(x_axis)
        pdict['y'].append(y_axis)
        pdict['z'].append(z_axis)

    pdict['N'] = np.stack(pdict['N'])
    pdict['Ca'] = np.stack(pdict['Ca'])
    pdict['C'] = np.stack(pdict['C'])
    pdict['Cb'] = np.stack(pdict['Cb'])

    # local frame for protein residues
    z_protein = pdict['Cb'][:protein_len] - pdict['Ca'][:protein_len]
    z_protein /= np.linalg.norm(z_protein, axis=-1)[:,None]
    x_protein = np.cross(pdict['Ca'][:protein_len]-pdict['N'][:protein_len], z_protein)
    x_protein /= np.linalg.norm(x_protein, axis=-1)[:,None]
    y_protein = np.cross(z_protein, x_protein)
    y_protein /= np.linalg.norm(y_protein, axis=-1)[:,None]

    if hetatm_residues:
        # local frame for heteroatoms
        z_hetatm = np.stack(pdict['z'])
        x_hetatm = np.stack(pdict['x'])
        y_hetatm = np.stack(pdict['y'])

        # combining protein and heteroatom local frames
        z = np.concatenate([z_protein, z_hetatm])
        x = np.concatenate([x_protein, x_hetatm])
        y = np.concatenate([y_protein, y_hetatm])
    else:
        x,y,z = x_protein, y_protein, z_protein

    xyz = np.stack([x,y,z])

    pdict['lfr'] = np.transpose(xyz, [1,0,2])

    start, end, j = 0, 0, 0
    atom_idx = [-1 for _ in range(atom_xyz.shape[0])]
    for i in range(len(atom_lst)):
        start = end
        end += atom_lst[i]
        atom_idx[start:end] = [j]*atom_lst[i]
        j = j+1

    def process_atoms(i):
        res_idx = atom_idx[i]
        lfr = pdict['lfr'][res_idx]
        p_local = np.einsum('ij,kj->ki', lfr, atom_xyz - atom_xyz[i])

        q_local = np.einsum('ij,kj->ki', lfr, pdict['lfr'][atom_idx][:, 0, :])
        k_local = np.einsum('ij,kj->ki', lfr, pdict['lfr'][atom_idx][:, 1, :])
        t_local = np.einsum('ij,kj->ki', lfr, pdict['lfr'][atom_idx][:, 2, :])

        return p_local, q_local, k_local, t_local

    results = Parallel(n_jobs=1)(delayed(process_atoms)(i) for i in range(atom_xyz.shape[0]))

    # Aggregate results just after processing them in parallel.
    p, q, k, t = zip(*results)
    p = np.stack(p)
    q = np.stack(q)
    k = np.stack(k)
    t = np.stack(t)

    return p, q, k, t

