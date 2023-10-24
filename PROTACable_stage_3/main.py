"""
------------------------------------------------------------------------------------------
Copyright (C) 2023 Hazem Mslati. All Rights Reserved.

This script is provided "AS IS" without any warranty.
Unauthorized copying, modification, or distribution is prohibited.
------------------------------------------------------------------------------------------
"""

import os
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm
from elabmol import ElabMols
import pymol
from multiprocessing import Pool, cpu_count
from openbabel import pybel
ob = pybel.ob
from openbabel import openbabel


def parse_arguments():
    parser = argparse.ArgumentParser(description='Process pdb files.')
    parser.add_argument('POCKET', type=str, help='Path to the POCKET pdb file.')
    parser.add_argument('POI_LIGAND', type=str, help='Path to the POI LIGAND pdb file.')
    parser.add_argument('E3_LIGAND', type=str, help='Path to the E3 LIGAND pdb file.')
    parser.add_argument('--LINKER', type=str, default=None, help='Path to the LINKER pdb file.')
    return parser.parse_args()

def merge_with_pymol(pocket_pdb_path, molecule_pdb_path, output_path):
    pymol.finish_launching(['pymol', '-qc'])  # '-qc' stands for "quiet and no GUI"
    pymol.cmd.delete("all")

    pymol.cmd.load(pocket_pdb_path, 'pocket')
    pymol.cmd.load(molecule_pdb_path, 'molecule')

    pymol.cmd.save(output_path, 'all')  # 'all' indicates we are saving everything loaded

#    pymol.cmd.quit()

def convert_pdb_to_mol2(pdb_input_path, mol2_output_path=None):
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdb", "mol2")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, pdb_input_path)
    mol.AddHydrogens()

    if not mol2_output_path:
        mol2_output_path = f"{os.path.splitext(pdb_input_path)[0]}.mol2"

    obConversion.WriteFile(mol, mol2_output_path)
    return mol2_output_path


def process_linker(linker_full_path, pocket, poi_ligand, e3_ligand, output_dir):
    em = ElabMols()
    result_mol = em.elaborate_two_mols(pocket, poi_ligand, e3_ligand, frag_path=linker_full_path)

    if result_mol is None:
        print(f"Linker {os.path.basename(linker_full_path)} was rejected")
        return

    # Temporary file to store the result molecule in PDB format
    temp_output_path = os.path.join(output_dir, "temp_result.pdb")
    Chem.MolToPDBFile(result_mol, temp_output_path)

    output_name_base = f"{os.path.splitext(os.path.basename(poi_ligand))[0]}_{os.path.splitext(os.path.basename(linker_full_path))[0]}"
    output_path_pdb = os.path.join(output_dir, output_name_base + ".pdb")

    # Merging with pymol
    merge_with_pymol(pocket, temp_output_path, output_path_pdb)

    ligands_dir = os.path.join(output_dir, 'ligands')
    if not os.path.exists(ligands_dir):
        os.makedirs(ligands_dir)

    # Save the resulting molecule as mol2 and sdf inside the 'ligands' directory
    output_path_mol2 = os.path.join(ligands_dir, output_name_base + ".mol2")
    output_path_sdf = os.path.join(ligands_dir, output_name_base + ".sdf")

    # Use the conversion function to get mol2 from the temporary pdb
    convert_pdb_to_mol2(temp_output_path, output_path_mol2)
    AllChem.MolToMolFile(result_mol, output_path_sdf)

    os.remove(temp_output_path)
    print(f"Linker {os.path.basename(linker_full_path)} was accepted. Results saved in {output_dir}")

def main():
    args = parse_arguments()

    output_dir = "output_stage_3_4"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    print(f"Output directory created at: {os.path.abspath(output_dir)}")

    linkers_to_process = []
    if args.LINKER:
        linkers_to_process.append(args.LINKER)
    else:
        # Looping over the collected linkers
        for linker_path in os.listdir(f"{os.environ['PROTACable']}/PROTACable_stage_3/linker_library/pdbs/"):
            if linker_path.endswith('.pdb'):
                linkers_to_process.append(os.path.join(f"{os.environ['PROTACable']}/PROTACable_stage_3/linker_library/pdbs/", linker_path))

    # Processing linkers using a simple loop
    for linker in tqdm(linkers_to_process):
        process_linker(linker, args.POCKET, args.POI_LIGAND, args.E3_LIGAND, output_dir)

if __name__ == "__main__":
    main()

