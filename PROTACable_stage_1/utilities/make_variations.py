"""
------------------------------------------------------------------------------------------
Copyright (C) 2023 Hazem Mslati. All Rights Reserved.

This script is provided "AS IS" without any warranty.
Unauthorized copying, modification, or distribution is prohibited.
------------------------------------------------------------------------------------------
"""


from rdkit import Chem
from rdkit.Chem import AllChem
import sys
import argparse

def extend_with_group(mol, atom_idx, group_smiles):
    """
    Extend molecule at a specific atom index with a given group.
    """
    group = Chem.MolFromSmiles(group_smiles)
    AllChem.EmbedMolecule(group, AllChem.ETKDG())
    
    combo = Chem.CombineMols(mol, group)
    edcombo = Chem.EditableMol(combo)
    edcombo.AddBond(atom_idx, mol.GetNumAtoms(), order=Chem.rdchem.BondType.SINGLE)
    new_mol = edcombo.GetMol()

    Chem.SanitizeMol(new_mol)
    
    # Optimize only the newly added group without disturbing the rest
    conf = new_mol.GetConformer(0)
    ff = AllChem.MMFFGetMoleculeForceField(new_mol, AllChem.MMFFGetMoleculeProperties(new_mol), confId=0)
    
    # Fix every atom except the newly added ones
    for i in range(mol.GetNumAtoms()):
        ff.AddFixedPoint(i)

    try:
        ff.Minimize(maxIts=1000)
    except:
        print("Warning: Could not optimize molecule.")

    return new_mol

def main(args):
    mol = Chem.MolFromPDBFile(args.pdb_file)
    if not mol:
        print(f"Error: Couldn't read molecule from {args.pdb_file}")
        sys.exit(1)
    
    Chem.SanitizeMol(mol)

    original_out_name = f"{args.output_prefix}_original.pdb"
    Chem.MolToPDBFile(mol, original_out_name)

    carboxyl_mol = extend_with_group(mol, args.exit_vector - 1, "C(=O)O")
    carboxyl_out_name = f"{args.output_prefix}_carboxyl.pdb"
    Chem.MolToPDBFile(carboxyl_mol, carboxyl_out_name)

    amide_mol = extend_with_group(mol, args.exit_vector - 1, "C(=O)N")
    amide_out_name = f"{args.output_prefix}_amide.pdb"
    Chem.MolToPDBFile(amide_mol, amide_out_name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extend ligand with amide and carboxyl groups.")
    parser.add_argument("pdb_file", help="Input PDB file")
    parser.add_argument("exit_vector", type=int, help="Exit vector for extension")
    parser.add_argument("output_prefix", help="Output file prefix")

    args = parser.parse_args()

    main(args)

