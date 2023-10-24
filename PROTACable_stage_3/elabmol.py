"""
------------------------------------------------------------------------------------------
Copyright (C) 2023 Hazem Mslati. All Rights Reserved.

This script is provided "AS IS" without any warranty.
Unauthorized copying, modification, or distribution is prohibited.
------------------------------------------------------------------------------------------
"""

import os
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdmolops
from rdkit.Chem import rdMolTransforms as rdmt
import openbabel
import copy
from openbabel import pybel

ob = pybel.ob

import warnings
warnings.simplefilter(action='ignore', category=UserWarning)

from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

class ElabMols():

    def rdkitmol_to_obmol(self, rdkitmol):
        mol = Chem.MolToMolBlock(rdkitmol)
        obConversion = openbabel.OBConversion()
        obConversion.SetInFormat("mol")
        obmol = openbabel.OBMol()
        obConversion.ReadString(obmol, mol)
        return obmol

    def get_protein_clashes(self, pdb, conformers):
        # Load the PDB file using Pybel
        #pdb = next(pybel.readfile('pdb', pdb_file))
        protein = pybel.Molecule(pdb)

        # Get the atom coordinates for the protein
        protein_coords = np.array([atom.coords for atom in protein.atoms])

        # Get the atom coordinates for each conformer
        atom_coords = []
        for conformer in conformers:
            # Convert RDKit mol object to OBMol
            mol = self.rdkitmol_to_obmol(conformer)
            pybel_mol = pybel.Molecule(mol)

            coords = np.array([atom.coords for atom in pybel_mol.atoms])

            atom_coords.append(coords)

        # Check for clashes between all pairs of atoms
        clash_counts = []
        for i in range(len(conformers)):
            coords1 = atom_coords[i]
            distances = np.sqrt(np.sum((protein_coords[:, None, :] - coords1[None, :, :])**2, axis=-1))
            clashes = np.where(distances < 1.0)

            # Count the number of clashes
            clash_counts.append(len(clashes[0]))

        return sum(clash_counts)

    def get_optimal_conformer(self, pdb_file, conformers, return_clashes=False):
        # Load the PDB file using Pybel
        pdb = next(pybel.readfile('pdb', pdb_file))
        protein = pybel.Molecule(pdb)

        # Get the atom coordinates for the protein
        protein_coords = np.array([atom.coords for atom in protein.atoms])

        # Get the atom coordinates for each conformer
        atom_coords = []
        for conformer in conformers:
            # Convert RDKit mol object to OBMol
            mol = self.rdkitmol_to_obmol(conformer)
            pybel_mol = pybel.Molecule(mol)

            coords = np.array([atom.coords for atom in pybel_mol.atoms])

            atom_coords.append(coords)

        # Check for clashes between all pairs of atoms 
        clash_counts = []
        for i in range(len(conformers)):
            coords1 = atom_coords[i]
            distances = np.sqrt(np.sum((protein_coords[:, None, :] - coords1[None, :, :])**2, axis=-1))
            clashes = np.where(distances < 1.0)

            # Count the number of clashes 
            clash_counts.append(len(clashes[0]))

        # Find the index
        optimal_index = np.argmin(clash_counts)

        if return_clashes:
            return clash_counts[optimal_index]
        else:
            return optimal_index
    
    def connectMols(self, mol1, mol2, neigh1, neigh2, radiolabels):
        #AllChem.EmbedMolecule(mol2)
        combined = Chem.CombineMols(mol1, mol2)
        emol = Chem.EditableMol(combined)
        bond_order = mol2.GetAtomWithIdx(neigh2[0]).GetBonds()[0].GetBondType()

        # Initial connection
        atom1 = mol1.GetAtomWithIdx(neigh1)
        atom2 = mol2.GetAtomWithIdx(neigh2[1])
        neighbor1_idx = atom1.GetNeighbors()[0].GetIdx()
        neighbor2_idx = atom2.GetNeighbors()[0].GetIdx()
        emol.AddBond(neighbor1_idx, neighbor2_idx + mol1.GetNumAtoms(), order=bond_order)
        emol.RemoveAtom(atom2.GetIdx() + mol1.GetNumAtoms())
        emol.RemoveAtom(atom1.GetIdx())

        # Additional connection if length of neigh2 is 4
        if '92' in radiolabels and len(neigh2) > 4:
            atom1 = mol1.GetAtomWithIdx(neigh1)
            atom2 = mol2.GetAtomWithIdx(neigh2[3])  # Connecting the 3rd atom of du2
            neighbor1_idx = atom1.GetNeighbors()[0].GetIdx()
            neighbor2_idx = atom2.GetNeighbors()[0].GetIdx()
            emol.AddBond(neighbor1_idx, neighbor2_idx + mol1.GetNumAtoms(), order=bond_order)
            emol.RemoveAtom(atom2.GetIdx() + mol1.GetNumAtoms())
            emol.RemoveAtom(atom1.GetIdx())

        mol = emol.GetMol()
 #       att_94 = neigh94 + mol1.GetNumAtoms()
        return mol  #, att_94

    def connectMols_e3(self, mol1, mol2, neigh1, neigh2):
        combined = Chem.CombineMols(mol1, mol2)
        emol = Chem.EditableMol(combined)
        bond_order = mol2.GetAtomWithIdx(neigh2[0]).GetBonds()[0].GetBondType()

        atom1 = mol1.GetAtomWithIdx(neigh1[1])
        atom2 = mol2.GetAtomWithIdx(neigh2[1])
        neighbor1_idx = atom1.GetNeighbors()[0].GetIdx()
        neighbor2_idx = atom2.GetNeighbors()[0].GetIdx()
        emol.AddBond(neighbor1_idx, neighbor2_idx + mol1.GetNumAtoms(), order=bond_order)
        emol.RemoveAtom(atom2.GetIdx() + mol1.GetNumAtoms())
        emol.RemoveAtom(atom1.GetIdx())

        mol = emol.GetMol()
        return mol


    @staticmethod
    def replace_dummy_with_hydrogen(mol):
        dummy_indices = []
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == '*':
                dummy_indices.append(atom.GetIdx())
                atom.SetAtomicNum(1)  # Replace with hydrogen
        return dummy_indices

    @staticmethod
    def replace_hydrogen_with_dummy(mol, dummy_indices):
        for idx in dummy_indices:
            mol.GetAtomWithIdx(idx).SetAtomicNum(0)  # Replace back with dummy atom
            
    @staticmethod
    def adjust_bond_type(mols, atom_idx1, atom_idx2, bond_type):
        for mol in mols:
            bond = mol.GetBondBetweenAtoms(atom_idx1, atom_idx2)
            bond.SetBondType(bond_type)

    def align(self, protein, mol1, mol2, mol3, du1, du2, ca1, ca2, ca3, radiolabels):
        du2_ = du2[1]
        ca2_ = ca2[0]
 
        neigh94 = ca2[-1]

        # Setup embedding parameters
        param = rdDistGeom.ETKDGv2()
        param.pruneRmsThresh = 0.1

        cids = rdDistGeom.EmbedMultipleConfs(mol2, 1000, param)

        #writer = Chem.SDWriter('output.sdf')
        #for conf_id in mol2.GetConformers():
        #    tmp_mol = Chem.Mol(mol2, False, conf_id.GetId())  # Creates a shallow copy with only the specific conformer
        #    writer.write(tmp_mol)
        #writer.close()


        alpha = 1  
        beta = 1  

        best_score = float('inf')  # the lower the score, the better
        best_distance = float('inf')
        final_molecule = None

        for cid in cids:

            mol2_conf = copy.deepcopy(mol2)
            mol2_conf.RemoveAllConformers()
            mol2_conf.AddConformer(mol2.GetConformer(cid))

            molIndex1 = mol1.GetAtomWithIdx(du1).SetAtomicNum(1)
            molIndex2 = mol2_conf.GetAtomWithIdx(du2_).SetAtomicNum(1)

            aligned_mol = Chem.Mol(mol2_conf.ToBinary())
            aligned = Chem.rdMolAlign.AlignMol(aligned_mol, mol1, atomMap=((ca2[0], du1), (du2[1], ca1)))
            connected = self.connectMols(mol1, aligned_mol, du1, du2, radiolabels)

            protein_clashes = self.get_protein_clashes(protein, [aligned_mol])
            pybel_conformer = pybel.Molecule(self.rdkitmol_to_obmol(aligned_mol))
            intramolecular_clashes = self.get_intermolecular_clashes(pybel.Molecule(self.rdkitmol_to_obmol(mol1)), pybel_conformer)
            intra_clashes_mol3 = self.get_intermolecular_clashes(pybel.Molecule(self.rdkitmol_to_obmol(mol3)), pybel_conformer)
            total_clashes = protein_clashes + intramolecular_clashes + intra_clashes_mol3

            distance = np.linalg.norm(np.array(aligned_mol.GetConformer().GetAtomPosition(ca2[-2])) - np.array(mol3.GetConformer().GetAtomPosition(ca3)))

            score = alpha * total_clashes + beta * distance

            if score < best_score:
                best_score = score
                best_distance = distance
                print(f"Minimizing intermolecular and intramolecular clash whil also minimizing distance between Linker and E3 ligase ligand: {distance}")
                final_molecule = connected

        return final_molecule, best_distance

    def rotate_dihedral(self, mol, begin_atom_idx, end_atom_idx, angle):
        """
        Rotates part of a molecule around a bond to achieve a desired dihedral angle.
        """
        conf = mol.GetConformer()

        # Atom positions
        begin_atom_pos = conf.GetAtomPosition(begin_atom_idx)
        end_atom_pos = conf.GetAtomPosition(end_atom_idx)

        # Create rotation axis
        rotation_axis = np.array([end_atom_pos.x - begin_atom_pos.x, 
                                  end_atom_pos.y - begin_atom_pos.y, 
                                  end_atom_pos.z - begin_atom_pos.z])
        rotation_axis /= np.linalg.norm(rotation_axis)

        # Define atoms to be moved during rotation
        # Here, we use a simple BFS approach to select the part of the molecule to rotate.
        atoms_to_move = set()
        visited = set()
        queue = [end_atom_idx]
        visited.add(end_atom_idx)
        while queue:
            current_atom = queue.pop()
            atoms_to_move.add(current_atom)
            for neighbor in mol.GetAtomWithIdx(current_atom).GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited:
                    visited.add(neighbor_idx)
                    queue.append(neighbor_idx)

        # Compute rotation matrix
        theta = np.radians(angle)
        a = rotation_axis
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)
        rotation_matrix = np.array([
            [cos_theta + a[0]**2 * (1 - cos_theta), a[0] * a[1] * (1 - cos_theta) - a[2] * sin_theta, a[0] * a[2] * (1 - cos_theta) + a[1] * sin_theta],
            [a[1] * a[0] * (1 - cos_theta) + a[2] * sin_theta, cos_theta + a[1]**2 * (1 - cos_theta), a[1] * a[2] * (1 - cos_theta) - a[0] * sin_theta],
            [a[2] * a[0] * (1 - cos_theta) - a[1] * sin_theta, a[2] * a[1] * (1 - cos_theta) + a[0] * sin_theta, cos_theta + a[2]**2 * (1 - cos_theta)]
        ])

        # Apply the rotation to the specified atoms
        for atom_idx in atoms_to_move:
            pos = conf.GetAtomPosition(atom_idx)
            relative_pos = np.array([pos.x - begin_atom_pos.x, pos.y - begin_atom_pos.y, pos.z - begin_atom_pos.z])
            new_pos = rotation_matrix.dot(relative_pos)
            new_pos += np.array([begin_atom_pos.x, begin_atom_pos.y, begin_atom_pos.z])
            conf.SetAtomPosition(atom_idx, new_pos)



    def identify_rotatable_bonds(self, mol):
        rotatable_bond = Chem.MolFromSmarts("[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]")
        return mol.GetSubstructMatches(rotatable_bond)

    def generate_conformers_by_torsion(self, mol, granularity=10):
        orig_mol = Chem.Mol(mol)
        mol.RemoveAllConformers()
        mol = Chem.AddHs(mol)
        cids = [mol.AddConformer(Chem.Conformer(orig_mol.GetConformer()), assignId=True)]

        rotatable_bonds = self.identify_rotatable_bonds(mol)
        for bond_indices in rotatable_bonds:
            print(bond_indices)
            begin_angle = 0
            end_angle = 360
            num_steps = int((end_angle - begin_angle) / granularity)
            for step in range(num_steps):
                angle = begin_angle + step * granularity
                new_conf = Chem.Conformer(mol.GetConformer(cids[0]))
                rdmt.SetDihedralDeg(new_conf, bond_indices[0], bond_indices[1], bond_indices[2], bond_indices[3], angle)
                mol.AddConformer(new_conf, assignId=True)

        return mol 

    @staticmethod
    def pdb_to_mol(pdb_path, use_rdkit=False):
        if use_rdkit:
            return Chem.MolFromPDBFile(pdb_path)
        else:
            mol = next(pybel.readfile('pdb', pdb_path))
            return mol.OBMol

    def get_attachment_indices_from_rdkit_mol(self, mol, attachment_points):
        # Sorting attachment points to ensure order
        attachment_points = sorted(attachment_points, key=lambda x: (int(x[-2:]), x[:-2]))

        att_indices = []
        dbs_idx = []
        radiolabels = []

        # Used to keep track of which attachment points we've processed
        processed_points = set()

        for point in attachment_points:
            for atom in mol.GetAtoms():
                if atom.GetPDBResidueInfo():
                    atom_name = atom.GetPDBResidueInfo().GetName().strip()
                    if atom_name == point:
                        if atom_name[-2:] not in processed_points:
                            radiolabels.append(atom_name[-2:])
                            att_indices.append(atom.GetIdx())

                            # Check for neighbor hydrogens
                            for nei in atom.GetNeighbors():
                                if nei.GetAtomicNum() == 1:
                                    att_indices.append(nei.GetIdx())

                            processed_points.add(atom_name[-2:])

                        # Check for the dbs_idx condition
                        if int(atom_name[-2:]) == 93:
                            dbs_idx.append(atom.GetIdx())

        return att_indices, dbs_idx, radiolabels

    def replace_with_dummy(self, mol, indices):
        for idx in indices:
            mol.GetAtomWithIdx(idx).SetAtomicNum(0)  # Replace hydrogen with dummy atom
            
    def get_intermolecular_clashes(self, mol1, mol2):
        mol1_coords = np.array([atom.coords for atom in mol1.atoms])
        mol2_coords = np.array([atom.coords for atom in mol2.atoms])

        distances = np.sqrt(np.sum((mol1_coords[:, None, :] - mol2_coords[None, :, :])**2, axis=-1))
        clashes = np.where(distances < 1.0)

        return len(clashes[0])

    def rotate_molecule(self, mol, theta, cid):
        # Define a rotation matrix for 60 degrees about the z-axis
        rotation_matrix = np.array([
            [np.cos(theta), -np.sin(theta), 0],
            [np.sin(theta), np.cos(theta),  0],
            [0,             0,              1]
        ])

        conf = mol.GetConformer(cid)
        for atom_idx in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(atom_idx)
            new_pos = rotation_matrix.dot(np.array([pos.x, pos.y, pos.z]))
            conf.SetAtomPosition(atom_idx, new_pos)
        return mol

    def elaborate_two_mols(self, protein_path, mol1_path, mol2_path, frag_path=None, double_bond_flags=None):
        protein = self.pdb_to_mol(protein_path)
        mol1 = self.pdb_to_mol(mol1_path, use_rdkit=True)
        mol1=AllChem.AddHs(mol1,addCoords=True)
        mol2 = self.pdb_to_mol(mol2_path, use_rdkit=True)
        mol2=AllChem.AddHs(mol2,addCoords=True)

        mol1_attachment_points = ["C90", "N90", "O90", "S90", "P90", "C94", "N94", "O94", "S94", "P94"]
        mol1_indices, _ , _ = self.get_attachment_indices_from_rdkit_mol(mol1, mol1_attachment_points)

        mol2_attachment_points = ["C95", "N95", "O95", "S95", "P95"]
        mol2_indices, _ , _ = self.get_attachment_indices_from_rdkit_mol(mol2, mol2_attachment_points)

        final_molecule = None
        if frag_path:  # If a fragment is provided
            frag = self.pdb_to_mol(frag_path, use_rdkit=True)
            frag = AllChem.AddHs(frag,addCoords=True)
            frag_attachment_points = ["C91", "N91", "O91", "S91", "P91",
                                      "C92", "N92", "O92", "S92", "P92",
                                      "C93", "N93", "O93", "S93", "P93",
                                      "C94", "N94", "O94", "S94", "P94"]
            frag_indices, dbs_idx, radiolabels = self.get_attachment_indices_from_rdkit_mol(frag, frag_attachment_points)
            print("radiolabels: ", radiolabels)
            double_bond_flags = dbs_idx

        #    aligned_frag_to_mol1_conformers = self.align(mol1, frag, mol1_indices[1], frag_indices, mol1_indices[0], frag_indices)
            final_molecule, best_distance = self.align(protein, mol1, frag, mol2, mol1_indices[1], frag_indices, mol1_indices[0], frag_indices, mol2_indices[0], radiolabels)
            discard = False
            if best_distance > 20.0:
                discard = True
                final_molecule = None
                return final_molecule

            if double_bond_flags and len(frag_indices) == 4:
                self.adjust_bond_type([final_molecule], mol1_indices[0], frag_indices[2], Chem.BondType.DOUBLE)
            elif double_bond_flags and len(frag_indices) == 6:
                self.adjust_bond_type([final_molecule], mol1_indices[0], frag_indices[4], Chem.BondType.DOUBLE)
            elif double_bond_flags and len(frag_indices) == 8:
                self.adjust_bond_type([final_molecule], mol1_indices[0], frag_indices[6], Chem.BondType.DOUBLE)
           
            e3_attachment_points = ["C94", "N94", "O94", "S94", "P94"]
            att94_indices, _ , _ = self.get_attachment_indices_from_rdkit_mol(final_molecule, e3_attachment_points)

            protac = self.connectMols_e3(final_molecule, mol2, att94_indices, mol2_indices)
            #elif double_bond_flags:
            #    self.adjust_bond_type([final_molecule], mol1_indices[0], frag_indices[0], Chem.BondType.DOUBLE)

        else:  # If no fragment and only mol1 has attachment for mol2
            pass

        return protac 


