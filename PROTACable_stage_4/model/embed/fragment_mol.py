"""
------------------------------------------------------------------------------------------
Copyright (C) 2023 Hazem Mslati. All Rights Reserved.

This script is provided "AS IS" without any warranty.
Unauthorized copying, modification, or distribution is prohibited.
------------------------------------------------------------------------------------------
"""

# fragment_mol.py

import os
from rdkit import Chem
from rdkit.Chem import AllChem, HybridizationType
import numpy as np
import itertools
import string
import copy
from Bio.PDB import PDBIO, PDBParser, Select
import tempfile
import random

class ResidueRenumber(Select):
    def __init__(self, new_resnames):
        self.new_resnames = new_resnames

    def accept_residue(self, residue):
        # Set new residue name
        residue.resname = self.new_resnames[residue.id[1]]
        return 1  # Accept the entity


class FragmentPDB:
    def __init__(self, protacs_embedding=False):
        self.protacs_embedding=protacs_embedding
        self.MOL_SPLIT_START = 57
        self.generated_resnames = set()
        self.hybridization_types = {HybridizationType.SP2, HybridizationType.SP}

    def unique_resnames(self):
        while True:
            resname = ''.join(random.choice(string.ascii_uppercase) for _ in range(3))
            if resname not in self.generated_resnames:
                self.generated_resnames.add(resname)
                return resname

    def okToBreak(self, bond):
        """
        Here we apply a bunch of rules to judge if the bond is OK to break.

        Parameters
        ----------
        bond : RDkit Bond
            The bond to consider breaking

        Returns
        -------
        Boolean : 
            OK or not to break.
        """

        # List of forbidden atoms
        forbidden_atoms = [15, 16]  # atomic numbers for phosphorus (P) and sulfur (S)

        # Do not break bonds that are part of a ring
        if bond.IsInRing():
            return False

        # Only break single bonds
        if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
            return False

        # Get the atoms at the beginning and end of the bond
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()

#         # Do not break bonds where either atom is in a ring
#         if begin_atom.IsInRing() or end_atom.IsInRing():
#             return False
        
        # this is kinda redundant but meh. Meh? You know what i mean..
        begin_atom_in_ring = begin_atom.IsInRing()
        end_atom_in_ring = end_atom.IsInRing()
        begin_atom_hybridization = begin_atom.GetHybridization()
        end_atom_hybridization = end_atom.GetHybridization()

        if (begin_atom_in_ring or end_atom_in_ring) and \
           (begin_atom_hybridization in self.hybridization_types and 
            end_atom_hybridization in self.hybridization_types):
            return False

        # Do not break bonds with hydrogen atoms
        if begin_atom.GetAtomicNum() == 1 or end_atom.GetAtomicNum() == 1:
            return False

        # Do not break bonds where either atom has atomic number >= MOL_SPLIT_START
        if begin_atom.GetAtomicNum() >= self.MOL_SPLIT_START or end_atom.GetAtomicNum() >= self.MOL_SPLIT_START:
            return False
        if self.protacs_embedding:

             # Do not break bonds with forbidden atoms
             if begin_atom.GetAtomicNum() in forbidden_atoms or end_atom.GetAtomicNum() in forbidden_atoms:
                 return False

#         # Do not break bonds where either atom is sp2 or sp hybridized
#         if begin_atom.GetHybridization() in [HybridizationType.SP2, HybridizationType.SP] or \
#            end_atom.GetHybridization() in [HybridizationType.SP2, HybridizationType.SP]:
#             return False

             # Do not break bonds where either atom is a terminal atom (atom with only one bond)
             if len(begin_atom.GetBonds()) == 1 or len(end_atom.GetBonds()) == 1:
                 return False
        
            # Do not break bonds where either atom has less than three heavy atom neighbors (for PCA)
             if len([neighbor for neighbor in begin_atom.GetNeighbors() if neighbor.GetAtomicNum() > 1]) < 4 or \
                len([neighbor for neighbor in end_atom.GetNeighbors() if neighbor.GetAtomicNum() > 1]) < 4:
                 return False

        elif not self.protacs_embedding:
            # Do not break bonds where either atom has less than three heavy atom neighbors (for PCA)
            if len([neighbor for neighbor in begin_atom.GetNeighbors() if neighbor.GetAtomicNum() > 1]) < 2 or \
               len([neighbor for neighbor in end_atom.GetNeighbors() if neighbor.GetAtomicNum() > 1]) < 2:
                return False


        # If none of the above rules apply, the bond can be broken
        return True

    def undo_id_label (self, frag, split_id):
        # I am trying to restore Hydrogens where the break happened
        for i, atom in enumerate(frag.GetAtoms()):
            if atom.GetAtomicNum() >= split_id:
                atom.SetAtomicNum(1)

        return frag

    def replace_dummy_with_hydrogen(self, mol):
        # Find dummy atoms
        dummy_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0]

        # Replace dummy atoms with hydrogen
        for i in dummy_indices:
            mol.GetAtomWithIdx(i).SetAtomicNum(1)  # Set atomic number to 1 for hydrogen
#         updated_mol = Chem.RemoveHs(mol)
        return mol

    def FragID_assign(self, mol):
        invariantID=AllChem.GetHashedMorganFingerprint(mol,radius=2,nBits=1024)
        key=str(''.join([str(item) for item in invariantID]))
        try:
            return FragID[key]
        except:
            FragID[key] = len(FragID)+1
            return [FragID[key]]
        
    def FragID_assign_2(self, mol):
            return list(MACCSkeys.GenMACCSKeys(mol))

    def is_empty(self, lst):
        if isinstance(lst, list):  # if it's a list
            return all((isinstance(sublst, list) and not sublst) or not isinstance(sublst, list) for sublst in lst) or not lst
        return False  # it's not a list

    # remove the fragment from the main molecule
    def remove_fragment(self, mol, frag):
        coords_frag=[]
        
        indices=[atom.GetIdx() for atom in frag.GetAtoms()]
        frag_conf=frag.GetConformer()
        for idx in indices:
            coords_frag.append(list(frag_conf.GetAtomPosition(idx)))

        to_remove=[]
        indices=[atom.GetIdx() for atom in mol.GetAtoms()]
        mol_conf=mol.GetConformer()
        for idx in indices: 
            if any(coord == tuple(round(y) for y in list(mol_conf.GetAtomPosition(idx))) \
                   for coord in [tuple(round(y) for y in x) for x in coords_frag]):
                to_remove.append(idx)

        # now remove the atoms on identified indices of fragments

        m = copy.deepcopy(mol)
        em1 = Chem.EditableMol(m)
        atomsToRemove=to_remove
        atomsToRemove.sort(reverse=True)
        for atom in atomsToRemove:
            em1.RemoveAtom(atom)

        m2 = em1.GetMol()
        Chem.SanitizeMol(m2)
        Chem.MolToSmiles(m2)

        return Chem.rdmolops.GetMolFrags(m2, asMols=True)

    def decomposed_atoms_idx(self, mol, coords_att):
        
        mol_idx=[]
        indices=[atom.GetIdx() for atom in mol.GetAtoms()]
        mol_conf=mol.GetConformer()
        for idx in indices: 
            if any(coord == tuple(round(y) for y in list(mol_conf.GetAtomPosition(idx))) \
                   for coord in [tuple(round(y) for y in x) for x in coords_att]):
                mol_idx.append(idx)
        return mol_idx
        
    # for hydrogen att points ---
    def hydrogen_frag_points(self, mol):
    
        hydrogenated = Chem.AddHs(mol)
        # find substructure hydrogen
        hyd_idx = hydrogenated.GetSubstructMatches(Chem.MolFromSmarts("[#1]"))
        # hydrogen indices are shifted by one
        hyd_idx_lst = list(x[0]+1 for x in hyd_idx) 

        ok_to_break = list()
        for i, atom in enumerate(hydrogenated.GetAtoms()):
            if i in hyd_idx_lst:
                ok_to_break.append(atom.GetBonds()[0])
        return [x.GetIdx() for x in ok_to_break]
    
    def get_atom_indices(self, original, fragment):
        """
        Retrieve atom indices in the original molecule based on the atom's coordinates in the fragment.

        Parameters
        ----------
        original : RDKit Mol object
            The original molecule before fragmentation.

        fragment : RDKit Mol object
            The molecule fragment.

        Returns
        -------
        atom_indices : List
            The atom indices in the original molecule for each atom in the fragment.
        """
        atom_indices = []
        orig_conf = original.GetConformer()
        frag_conf = fragment.GetConformer()

        for frag_idx in range(fragment.GetNumAtoms()):
            frag_pos = frag_conf.GetAtomPosition(frag_idx)

            for orig_idx in range(original.GetNumAtoms()):
                orig_pos = orig_conf.GetAtomPosition(orig_idx)

                if frag_pos.x == orig_pos.x and frag_pos.y == orig_pos.y and frag_pos.z == orig_pos.z:
                    atom_indices.append(orig_idx)
                    break

        return atom_indices
    

    def reorder_fragments_and_bonds(self, fragments, broken_bonds, mol):
        # Create a nested dictionary mapping from atoms to fragments using atomic coordinates
        atom_to_fragment = {}
        for fragment_idx, fragment in enumerate(fragments):
            conformer = fragment.GetConformer()
            for atom in fragment.GetAtoms():
                if atom.GetSymbol() == "*" : continue
                atom_index = atom.GetIdx()
                pos = conformer.GetAtomPosition(atom_index)
                coords = (pos.x, pos.y, pos.z)  
                atom_to_fragment[coords] = fragment_idx 

        # Define a function to get the fragment index based on a distance threshold
        def get_fragment_index(target_coords):
            for coords, fragment_idx in atom_to_fragment.items():
                if np.linalg.norm(np.array(coords) - np.array(target_coords)) < 1e-3:
                    return fragment_idx
            return None
        
        conformer = mol.GetConformer() 
        # Create a new list of bonds, where each bond is associated with two fragments
        fragment_bonds = []
        for bond in broken_bonds:
            atom1 = mol.GetAtomWithIdx(bond[0])
            atom2 = mol.GetAtomWithIdx(bond[1])
            pos1 = conformer.GetAtomPosition(atom1.GetIdx())
            pos2 = conformer.GetAtomPosition(atom2.GetIdx())
            coords1 = (pos1.x, pos1.y, pos1.z)  # changed to tuple from np.array
            coords2 = (pos2.x, pos2.y, pos2.z)  # changed to tuple from np.array
            fragment1 = get_fragment_index(coords1)
            fragment2 = get_fragment_index(coords2)
            if fragment1 is not None and fragment2 is not None and fragment1 != fragment2:
                fragment_bonds.append((fragment1, fragment2, bond))

#         # Sort the bonds by the fragments they connect
#         fragment_bonds.sort()

#         # Create new lists of fragments and broken_bonds in the correct order
#         new_fragments = [None] * len(fragments)
#         new_broken_bonds = [None] * len(broken_bonds)
#         for i, (fragment1, fragment2, bond) in enumerate(fragment_bonds):
#             new_fragments[i] = fragments[fragment1]
#             new_broken_bonds[i] = bond
#         new_fragments[-1] = fragments[fragment2]
#         # HERE WE MAY HAVE DUPLICATE FRAGS.
        return fragment_bonds


    # Divide a molecule into fragments
    def split_molecule(self, mols):
        accepted_fragments = {}
        original_mols = {}
        fragment_counts = []

        total_mols = 0
        for mol in mols:
            current_fragments = {}  # Initialize new dictionary for this molecule
            if not mol: continue
            split_id = self.MOL_SPLIT_START
            mol = Chem.AddHs(mol)
            total_mols+=1
            res = []
            ok_bonds_list = []

            to_check = [mol]
            while len(to_check) > 0:
                ok_bonds, ms = self.spf(mol, to_check.pop(), split_id)
                ok_bonds_list.append(ok_bonds)
                if len(ms) == 1:
                    res += ms
                else:
                    to_check += ms
                    split_id += 1

            if not self.is_empty(ok_bonds_list):
                ok_bonds = [number for group in ok_bonds_list for number in group]

                fragmented = Chem.FragmentOnBonds(mol, ok_bonds)
                fragments = list(Chem.GetMolFrags(fragmented, asMols=True))

                broken_bonds = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in mol.GetBonds() if bond.GetIdx() in ok_bonds]
                fragment_bonds = self.reorder_fragments_and_bonds(fragments, broken_bonds, mol)

                for fragment in fragments:
                    identifier = self.unique_resnames()
                    fragment = self.replace_dummy_with_hydrogen(fragment)
                    accepted_fragments[identifier] = fragment
                    current_fragments[identifier] = fragment

                fragment_counts.append([current_fragments, fragment_bonds, len(fragments)])
            else:
                identifier = self.unique_resnames()
                mol = self.replace_dummy_with_hydrogen(mol)
                accepted_fragments[identifier] = mol
                current_fragments[identifier] = mol

                fragment_counts.append([current_fragments, [], (1)])
           
            #original mol
            identifier = self.unique_resnames()
            mol = self.replace_dummy_with_hydrogen(mol)
            original_mols[identifier] = mol
            

        return accepted_fragments, fragment_counts, original_mols, total_mols

    # Function for doing all the nitty gritty splitting work.
    # loops over bonds until bonds get exhausted or bonds are ok to break, whichever comes first. If ok to break, then each
    # fragment needs to be checked individually again through the loop
    def spf(self, original, mol, split_id):
        
        ok_bonds=[]
        bonds = mol.GetBonds()
        for i in range(len(bonds)):
            if self.okToBreak(bonds[i]):
                
                mol = Chem.FragmentOnBonds(mol, [i])
                # Dummy atoms are always added last
                n_at = mol.GetNumAtoms()
#                 print('Split ID', split_id)
                att1=mol.GetAtomWithIdx(n_at-1)
                att2=mol.GetAtomWithIdx(n_at-2)
                
                mol_conf=mol.GetConformer()
                att1_idx=self.decomposed_atoms_idx(original, [list(mol_conf.GetAtomPosition(n_at-1))])
                att2_idx=self.decomposed_atoms_idx(original, [list(mol_conf.GetAtomPosition(n_at-2))])
                
                bond_idx=original.GetBondBetweenAtoms(att1_idx[0],att2_idx[0]).GetIdx()
                ok_bonds.append(bond_idx)
                
                att1.SetAtomicNum(split_id)
                att2.SetAtomicNum(split_id)
                return ok_bonds, Chem.rdmolops.GetMolFrags(mol, asMols=True)

        # If the molecule could not been split, return original molecule
        return [], [mol]

    def pdb_2_sdf(self, pdb):
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("pdb", "sdf")
        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, pdb)   # Open Babel will uncompress automatically

        mol.AddHydrogens()


        obConversion.WriteFile(mol, f"{pdb.split('.')[0]}.sdf")
        return f"{pdb.split('.')[0]}.sdf"
    
    def pdb_2_mol2(self, pdb):
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("pdb", "mol2")
        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, pdb)   # Open Babel will uncompress automatically

        mol.AddHydrogens()


        obConversion.WriteFile(mol, f"{pdb.split('.')[0]}.mol2")
        return f"{pdb.split('.')[0]}.mol2"
    
    def sdf_2_pdb(self, sdf):
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("sdf", "pdb")
        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, sdf)   # Open Babel will uncompress automatically

        mol.AddHydrogens()
        obConversion.WriteFile(mol, f"{sdf.split('.')[0]}.pdb")
        return f"{sdf.split('.')[0]}.pdb"

    def save_bpdb(self, pdb,ppdb, record):  
        ppdb.to_pdb(path=f"{record}_{pdb.split('.')[0].split('_')[0]}.pdb",
                    records=[record],
                    gz=False, 
                    append_newline=True)

#     def get_HOH_pdb(self, pdb):
#         ppdb = PandasPdb() 
#         ppdb.read_pdb(pdb) 
#         ppdb.df['HETATM']=ppdb.df['HETATM'].loc[ppdb.df['HETATM']['residue_name'].isin(self.values)]
#         ppdb.to_pdb(path=f"HOH_{pdb.split('.')[0].split('_')[0]}.pdb",
#                 records=['HETATM'],
#                 gz=False, 
#                 append_newline=True)

    def keep_relevant_hetatm(self, pdb):
        raw=str(self.pdb).split('/')[-1]
        with open(self.orig_pdb_drc, 'r') as f1, open(f"{raw.split('.')[0]}_protein.pdb", 'w') as f2:
            for line in f1:
                if 'ATOM' in line:
                    f2.write(line)
        with open(self.orig_pdb_drc, 'r') as f1, open(f"{raw.split('.')[0]}_ligand.pdb", 'w') as f2:
            for line in f1:
                if line.startswith('HETATM') and not any(ion in line for ion in self.ions):
                    putative_ligand=line[17:20]
                    print("putative ligand found: ", putative_ligand)
                    break
            for line in f1:
                if line.startswith('HETATM') and (putative_ligand in line) and not any(ion in line for ion in self.ions):
                    f2.write(line)
        return
    
    
    def fragment_and_plif(self):
        path = os.getcwd()
#         if not os.path.exists('fragment'):
#             os.mkdir(f'{path}/fragment')
#         os.chdir("fragment")

        raw=self.pdb.split('/')[-1].split('.')[0]
        self.keep_relevant_hetatm(self.pdb)
        lig_sdf=self.pdb_2_sdf(f'{raw}_ligand.pdb')
        lig_mol2=self.pdb_2_mol2(f'{raw}_ligand.pdb')
        
#         # in case the original pdb file is corrupted
#         content = open(f'{raw}_ligand.pdb').read()
#         hets=re.findall("^HETATM (.*)", content, re.M)
#         if len(hets)<5:

#             with open(f'{raw}_ligand.pdb', 'r') as file :
#                 filedata = file.read()

#             # Replace the target string
#             filedata = filedata.replace('ATOM  ', 'HETATM')

#             # Write the file out again
#             with open(f'{raw}_ligand.pdb', 'w') as file:
#                 file.write(filedata)
        
        # an extensive and elaborate multi-try methods to read into RdKit
        try: 
            fragment_mols = Chem.RemoveHs(fragment_mols[0])
            output_frags = self.split_molecule(fragment_mols,raw)
            
        except:  
            try:
                fragment_mols = Chem.SDMolSupplier(lig_sdf, removeHs=True, sanitize=False)
                output_frags = self.split_molecule(fragment_mols[0],raw)
            except:
                try: 
                    fragment_mols_alt = Chem.MolFromMol2File(lig_mol2, sanitize=True, removeHs=True)
                    output_frags = self.split_molecule(fragment_mols_alt,raw)
                except:
                    try: 
                        fragment_mols = AllChem.MolFromPDBFile(f'{raw}_ligand.pdb')
                        output_frags = self.split_molecule(fragment_mols,raw)
                    except:
                        raise Exception(f"Something wrong with ligand")
#         os.chdir(f'{path}')
        
        return output_frags
        
    def write_pdb(self, accepted_fragments, pdb_out_path):
        if os.path.exists(pdb_out_path):
            os.remove(pdb_out_path)
        for i, (identifier, frag) in enumerate(accepted_fragments.items()):
            pdb_block = Chem.MolToPDBBlock(frag)

            # Edit the residue name and exclude certain lines
            pdb_lines = []
            for line in pdb_block.split('\n'):
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    pdb_lines.append(line[:17] + '{:3s}'.format(identifier) + line[20:])
                    
                elif line.startswith('ENDHETATM'):
                    pdb_lines.append(line[3:])
                    
                elif not (line.startswith('CONECT') or line.startswith('END')):
                    pdb_lines.append(line)

            pdb_block = '\n'.join(pdb_lines)

            # Write the modified PDB block
            with open(pdb_out_path, 'a') as f:
                f.write(pdb_block)  # remove '\n' here

       
    def split_and_write(self, mols, pdb_out_path):
        accepted_fragments, fragment_counts, original_mols, total_mols = self.split_molecule(mols)
        self.write_pdb(accepted_fragments, f"{pdb_out_path}.pdb")
        self.write_pdb(original_mols, f"{pdb_out_path}_original_mols.pdb")
        return accepted_fragments, fragment_counts, total_mols


