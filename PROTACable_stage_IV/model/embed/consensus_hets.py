"""
------------------------------------------------------------------------------------------
Copyright (C) 2023 Hazem Mslati. All Rights Reserved.

This script is provided "AS IS" without any warranty.
Unauthorized copying, modification, or distribution is prohibited.
------------------------------------------------------------------------------------------
"""

import pymol
from pymol import cmd
import argparse
import os

# Initialize PyMOL
pymol.finish_launching(['pymol', '-qc'])  # Quiet and no GUI

common_ions = ['ZN', 'FE', 'MG', 'MN', 'K', 'NA', 'CA']  # Common ions list

def get_ions_and_ligand(pdb_file, obj_name):
    cmd.load(pdb_file, obj_name)

    # Select common ions and give them unique names
    ion_names = []
    for i, ion in enumerate(common_ions):
        ion_name = f'{obj_name}_ion_{i}'
        # Check if ion exists in the molecule before creating a selection
        if cmd.count_atoms(f'{obj_name} and resn {ion}'):
            cmd.select(ion_name, f'{obj_name} and resn {ion}')
            ion_names.append(ion_name)

    # Select largest HETATM ligand
    ligands = cmd.get_names('objects', 1, f'({obj_name} and hetatm)')
    largest_ligand = max(ligands, key=lambda name: cmd.count_atoms(f'{name} and hetatm'))

    # Rename largest ligand to a generic name
    cmd.set_name(largest_ligand, f'{obj_name}_largest_ligand')

    return f'{obj_name}_largest_ligand', ion_names

def main():
    parser = argparse.ArgumentParser(description='Process PDB files to keep only the largest ligand and common ions.')
    parser.add_argument('pdb_file', type=str, help='PDB file to be processed.')
    parser.add_argument('true_pdb', type=str, help='True PDB file.')

    args = parser.parse_args()

    # Get ions and largest ligand from each PDB
    ligand1, ions1 = get_ions_and_ligand(args.pdb_file, 'molecule1')
    ligand2, ions2 = get_ions_and_ligand(args.true_pdb, 'molecule2')

    # Determine common ions
    common_ion_selections = list(set(ions1) & set(ions2))

    # Save selections to new PDB files for each input
    for pdb_file, ligand in zip([args.pdb_file, args.true_pdb], [ligand1, ligand2]):
        selections = ' or '.join(common_ion_selections + [ligand])
        cmd.select('output', selections)
        cmd.save(pdb_file.replace('.pdb', '_filtered.pdb'), 'output')

        # Delete selections related to current molecule
        cmd.delete('sele*')

if __name__ == "__main__":
    main()

