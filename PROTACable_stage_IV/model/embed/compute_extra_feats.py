"""
------------------------------------------------------------------------------------------
Copyright (C) 2023 Hazem Mslati. All Rights Reserved.

This script is provided "AS IS" without any warranty.
Unauthorized copying, modification, or distribution is prohibited.
------------------------------------------------------------------------------------------
"""

import os
import torch
import esm
import joblib
import glob
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import three_to_one, is_aa
from tqdm import tqdm
import os
import joblib
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from sklearn.preprocessing import StandardScaler
from descriptastorus.descriptors.DescriptorGenerator import MakeGenerator
from openbabel import openbabel
from openbabel import pybel
import pymol
from sklearn.preprocessing import StandardScaler

# Define your scaler here (with .fit method)
#scaler = StandardScaler()


# Load ESM-2 model
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()  # disables dropout for deterministic results

# Parse protein sequences from PDB files
parser = PDBParser(QUIET=True)  # QUIET=True suppresses warnings



substitutions = {
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
    'TPL':'TRP', 'TPO':'THR', 'TPQ':'ALA', 'TRG':'LYS', 'TRO':'TRP', 'TYB':'TYR', 'TYI':'TYR', 'TYQ':'TYR', 'TYS':'TYR', 'TYY':'TYR'
}


def get_sequence_from_pdb_file_4_esm(file_path):
    parser = PDBParser()
    structure = parser.get_structure('protein', file_path)
    seq = ""
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname()
                if resname in substitutions:
                    resname = substitutions[resname]
                if not is_aa(resname, standard=True):
                    resname = 'GLY'  # Replace with Glycine if it's not a standard residue
                seq += three_to_one(resname)
    return seq

def compute_esm_embedding(file_path, batch_converter, model, alphabet):
    try:
        protein_sequence = get_sequence_from_pdb_file_4_esm(file_path)
    except Exception as e:
        print(f"Failed to parse {file_path}. Reason: {str(e)}")
        return None

    file_name = os.path.basename(file_path)
    data = [(file_name[:-4], protein_sequence)]  # create a tuple with protein name and sequence
    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

    # Extract per-residue representations (on CPU)
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[33])
    token_representations = results["representations"][33]

    # Generate per-sequence representations via averaging
    sequence_representations = token_representations[0, 1 : batch_lens[0] - 1].mean(0)

    # return the numpy array
    embedding_np = sequence_representations.cpu().numpy()
    return embedding_np


def save_ligand_with_pymol(pdb_path):
    pymol.cmd.load(pdb_path, 'molecule')
    pymol.cmd.remove('solvent')  # remove water and other solvent molecules if present
    pymol.cmd.select('ligand', 'hetatm')  # select all heteroatoms
    ligand_path = pdb_path.replace('.pdb', '_ligand.pymol.pdb')
    pymol.cmd.save(ligand_path, 'ligand')
    pymol.cmd.delete('all')  # clear PyMOL's memory
    return ligand_path

def pdb_2_mol2(pdb):
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdb", "mol2")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, pdb)   # Open Babel will uncompress automatically

    mol.AddHydrogens()

    mol2_path = pdb.replace('.pdb', '.mol2')
    obConversion.WriteFile(mol, mol2_path)
    return mol2_path

def pdb_to_mol2_with_pybel(pdb_path):
    mol = next(pybel.readfile('pdb', pdb_path))
    mol.addh()
    mol2_path = pdb_path.replace('.pdb', '_pybel.mol2')
    mol.write('mol2', mol2_path, overwrite=True)
    return mol2_path

def process_pdb_file(pdb_path, ligands_path, by_id):
    try:
        # Write a temporary HETATM file
        hetatm_file_path = pdb_path.replace('.pdb', '_HETATM.pdb')
        with open(pdb_path, 'r') as inp, open(hetatm_file_path, 'w') as out:
            for line in inp:
                if line.startswith('HETATM'):
                    out.write(line)

        mol=None
        if ligands_path is not None and by_id is not None:
            print(f"Loading ligands from {ligands_path} and searching for files with pattern {by_id}")
            mol2_path = glob.glob(ligands_path + f'/{by_id}*mol2')[0]
            sdf_path = glob.glob(ligands_path + f'/{by_id}*sdf')[0]
            if mol2_path and mol is None:
                mol = Chem.rdmolfiles.MolFromMol2File(mol2_path)
            if sdf_path and mol is None:
                mol = next(Chem.SDMolSupplier(sdf_path))

        # Calculate the Morgan fingerprints
        elif mol is None or ligands_path is None or by_id is None:
            print("Original ligands path directory is invalid. Retrying to extract ligand manually")
            pdb_block = open(hetatm_file_path).read()
            mol = Chem.rdmolfiles.MolFromPDBBlock(pdb_block)
        if mol is None:
            mol2_path = pdb_2_mol2(pdb_path)
            mol = Chem.rdmolfiles.MolFromMol2File(mol2_path)
        if mol is None:  # If RDKit still fails, try with Pybel
            pybel_mol2_path = pdb_to_mol2_with_pybel(pdb_path)
            mol = Chem.rdmolfiles.MolFromMol2File(pybel_mol2_path)
        if mol is None:  # If RDKit still fails, try with PyMOL
            pymol_ligand_path = save_ligand_with_pymol(pdb_path)
            mol = Chem.rdmolfiles.MolFromPDBFile(pymol_ligand_path)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        
        # Create the RDKit2D descriptor generator and process the SMILES
        generator = MakeGenerator(("rdkit2dnormalized",))
        smiles = Chem.MolToSmiles(mol)
        descriptors = generator.process(smiles)
        
        # Combine the Morgan fingerprints and RDKit2D descriptors
        features = np.concatenate((fp, descriptors[1:]))  # descriptors[0] is True/False for success/failure
        
        # Delete the temporary HETATM file
        os.remove(hetatm_file_path)

        return features
    
    except Exception as e:
        print(f'Failed to calculate features for {pdb_path}. Error: {e}')
        return None

def pdb_to_features(pdb_path, scaler_path, ligands_path=None, by_id=None):
    #scaler = joblib.load(scaler_path)
    # Process pdb and get features
    by_id = os.path.basename(by_id).split('.')[0]
    features = process_pdb_file(pdb_path, ligands_path, by_id)

    # Check if features were generated
    if features is None:
        print(f'Failed to generate features for {pdb_path}.')
        return None

    # Use the provided StandardScaler to scale the features
    features_scaled = features #scaler.transform(features.reshape(1, -1))

    return features_scaled


