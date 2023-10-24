"""
------------------------------------------------------------------------------------------
Copyright (C) 2023 Hazem Mslati. All Rights Reserved.

This script is provided "AS IS" without any warranty.
Unauthorized copying, modification, or distribution is prohibited.
------------------------------------------------------------------------------------------
"""

from Bio.PDB import *
from Bio import PDB
from Bio.Seq import Seq
from Bio import pairwise2
import sys, os, shutil
import tempfile
from subprocess import Popen, PIPE
import pymol2
import pymol
from pymol import cmd
import glob
import copy
from Bio.PDB import PDBParser, Superimposer, PDBIO, Atom, PDBIO, Select
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
import numpy as np
from scipy.optimize import linear_sum_assignment
from plif_utilities import SUBSTITUTIONS, clean_pdb



three_to_one = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

def robust_three_to_one(residue):
    """
    Convert a three-letter residue code to a one-letter code, 
    with a fallback for non-standard residues.
    """
    
    # If the residue is non-standard, get its standard equivalent
    if residue in SUBSTITUTIONS:
        residue = SUBSTITUTIONS[residue]

    # Return the one-letter code, or a placeholder if not found
    return three_to_one.get(residue, 'X')

def overwrite_pdb_with_pymol(pdb_file):
    """
    Reads and overwrites a PDB file using PyMOL.
    
    Parameters:
        pdb_file (str): Path to the PDB file to be overwritten.
    """
    
    # Start PyMOL without GUI
    pymol.finish_launching(['pymol', '-cqi'])
    
    # Load the PDB file into PyMOL
    pymol.cmd.load(pdb_file, "structure")
    
    # Save (overwrite) the PDB file from PyMOL
    pymol.cmd.save(pdb_file, "structure")
    

def uniquify_pdb(pdb_file):
    # PDB Parser
    parser = PDBParser(QUIET=True)

    chain_sequences = {}  # Dictionary to hold unique chains and their ligand and atom count

    # Load protein
    structure = parser.get_structure('protein', pdb_file)

    # New structure for unique chains
    new_structure = Structure.Structure('new_protein')

    # Create a new model for the structure
    model = Model.Model(0)

    # Polypeptide builder
    ppb = PPBuilder()

    # Iterate over each chain in the protein
    for chain1 in structure.get_chains():
        # Get the sequence for the chain
        if ppb.build_peptides(chain1):
            seq1 = ppb.build_peptides(chain1)[0].get_sequence()
        else:
            continue
        # Count of hetatms that are not water or ions
        ligand_atoms1 = [a for r in chain1.get_residues() if r.id[0].startswith("H") and r.resname not in ["HOH", "NA", "K", "MG", "CA", "ZN", "CL", "CD", "MN", "CO", "NI", "CU"] for a in r.get_atoms()]
#        #ligand_atoms1 = [a for a in chain1.get_atoms() if a.parent.id[0].startswith("H") and a.parent.resname not in ["HOH", "NA", "K", "MG", "CA", "ZN", "CL", "CD", "MN", "CO", "NI", "CU"]]
        #ligand_atoms1 = [a for r in chain1.get_residues() if r.id[0].startswith('H') for a in r.get_atoms()]

        ligand_count1 = len(ligand_atoms1)

        similar_chain_found = False
        for seq, (chain, ligand_count, ligand_atoms) in chain_sequences.items():
            # Align the sequences and calculate the similarity
            alignment = pairwise2.align.globalxx(str(seq1), str(seq), one_alignment_only=True, score_only=True)
            similarity = alignment / len(str(seq)) * 100

            # Check if the sequences are similar (90% similarity or higher)
            if similarity >= 90:
                similar_chain_found = True

                # If ligand count is equal, pick the chain with the highest total ligand atom count
                if ligand_count1 == ligand_count and len(ligand_atoms1) > len(ligand_atoms):
                    chain_sequences[seq] = (chain1, ligand_count1, ligand_atoms1)

                # If ligand count is higher than previous same sequence
                elif ligand_count1 > ligand_count:
                    chain_sequences[seq] = (chain1, ligand_count1, ligand_atoms1)

        # If chain sequence not seen before
        if not similar_chain_found:
            chain_sequences[seq1] = (chain1, ligand_count1, ligand_atoms1)

    # Add unique chains to the new model
    for seq, (chain, ligand_count, ligand_atoms) in chain_sequences.items():
        new_chain = Chain(chain.id)
        for residue in chain:
            new_chain.add(residue.copy())
        model.add(new_chain)

    # Add the model to the new structure
    new_structure.add(model)

    # Save the new structure with unique chains to a single PDB file
    io = PDBIO()
    io.set_structure(new_structure)
    io.save(os.path.splitext(pdb_file)[0] + '_unique_chains.pdb')
    return pdb_file.split('.')[0] + '_unique_chains.pdb'


# Define a selection class for Biopython to select only non-hetero atoms
class NonHetSelect(Select):
    def accept_atom(self, atom):
        return not atom.get_parent().get_resname() in ['HOH', 'WAT'] and not atom.is_disordered() != 0

# Start PyMOL
pm = pymol2.PyMOL()
pm.start()
cmd = pm.cmd

# Function to save a chain to a file
def save_chain(structure, chain_id, filename):
    class ChainSelect(Select):
        def accept_chain(self, chain):
            return chain.get_id() == chain_id

    io = PDBIO()
    io.set_structure(structure)
    io.save(filename, ChainSelect())

def remove_commands(file_name, ligand = False):
    with open(file_name, 'r') as file:
        lines = file.readlines()
    
    # Find the index of the first line that begins with "REMARK" or "ATOM"
    start_index = next(i for i, line in enumerate(lines) if line.startswith('REMARK') or line.startswith('ATOM'))

    # Keep only lines starting from the first "REMARK" or "ATOM"
    lines = lines[start_index:]
    filename = f"{file_name.split('.')[0]}_aligned_w_lig.pdb" if ligand else f"{file_name.split('.')[0]}_aligned.pdb"
    with open(filename, 'w') as file:
        file.writelines(lines)
        
def split_pdb_file(input_file, output_file1, output_file2):
    with open(input_file, 'r') as f, open(output_file1, 'w') as out1, open(output_file2, 'w') as out2:
        ter_found = False
        second_file = False
        for line in f:
            if ter_found and "ATOM" in line:
                second_file = True
                out2.write(line)
                continue
            if "TER" in line:
                ter_found = True
                continue
            if not second_file:
                out1.write(line)
            else:
                out2.write(line)

# Function to superimpose structures using TM-align
def superimpose(pdb1, pdb2):

    # Extract chains
    parser = PDBParser()
    structure1 = parser.get_structure("pdb1", pdb1)
    structure2 = parser.get_structure("pdb2", pdb2)

    for chain_obj1 in structure1.get_chains():
        for chain_obj2 in structure2.get_chains():
            
            # Create temporary files for each chain
            tmp1 = tempfile.mktemp('.pdb')
            tmp2 = tempfile.mktemp('.pdb')
            save_chain(structure1, chain_obj1.get_id(), tmp1)
            save_chain(structure2, chain_obj2.get_id(), tmp2)

            # Run TM-align
            output_filename = f'TM_{chain_obj1.get_id()}_{chain_obj2.get_id()}.sup'
            p = Popen(['TMalign', tmp1, tmp2, '-a', '-d', '8.0', '-o', output_filename], stdout=PIPE, stderr=PIPE)
            p.wait()

            # Read output
            out, err = p.communicate()
            lines = out.decode('utf-8').split('\n')

            # Only superimpose chains with sequence identity or overlap > 90%
            # Extract sequence identity
            seq_id_line = next(line for line in lines if 'Seq_ID' in line)
            seq_id = float(seq_id_line.split('=')[-1].strip())

            # Extract lengths of chains
            len_chain_1_line = next(line for line in lines if 'Length of Chain_1' in line)
            len_chain_1 = int(len_chain_1_line.split()[-2])

            len_chain_2_line = next(line for line in lines if 'Length of Chain_2' in line)
            len_chain_2 = int(len_chain_2_line.split()[-2])

            # Extract aligned length
            aligned_len_line = next(line for line in lines if 'Aligned length=' in line)
            aligned_len = int(aligned_len_line.split(',')[0].split('=')[-1].strip())

            # Calculate sequence overlap
            seq_overlap = max(aligned_len / len_chain_1, aligned_len / len_chain_2)

#             print(seq_id,seq_overlap )
            if seq_id >= 0.9 or seq_overlap >= 0.9:
                # Load superimposed structure into PyMOL
                if os.path.exists(f"TM_{chain_obj1.get_id()}_{chain_obj2.get_id()}.sup_all_atm_lig"):
                    remove_commands(f"TM_{chain_obj1.get_id()}_{chain_obj2.get_id()}.sup_all_atm_lig", ligand=True)
                    cleaned_pdb = f"TM_{chain_obj1.get_id()}_{chain_obj2.get_id()}_aligned_w_lig.pdb"
                else:
                    remove_commands(f"TM_{chain_obj1.get_id()}_{chain_obj2.get_id()}.sup_all_atm", ligand=False)
                    cleaned_pdb = f"TM_{chain_obj1.get_id()}_{chain_obj2.get_id()}_aligned_w_lig.pdb"
#                 print(chain1.get_id())
                # Now we need to separate the pdb files of tmalign output
                # Splitting the PDB files of tmalign output
                split_pdb_file(cleaned_pdb, f"TM_pdb1_{chain_obj1.get_id()}.pdb", f"TM_pdb2_{chain_obj2.get_id()}.pdb")
                
                cmd.load(f"TM_pdb1_{chain_obj1.get_id()}.pdb", 'sup1')
                chain_id1 = cmd.get_chains('sup1')[0]
                
                cmd.load(f"TM_pdb2_{chain_obj2.get_id()}.pdb", 'sup2')
                chain_id2 = cmd.get_chains('sup2')[0]
                
                cmd.align('sup1 and chain ' + str(chain_id1), 'sup2 and chain ' + str(chain_id2), object="aln")

    # Create the directory if it doesn't exist
    directory = "super"
    if not os.path.exists(directory):
        os.makedirs(directory)

    # Move the saved files to the "super" directory
    pdb1_filename = pdb1.split('/')[-1].split('.')[0]
    pdb2_filename = pdb2.split('/')[-1].split('.')[0]

    cmd.save(f"aligned_pymol_{pdb1_filename}.pdb", 'sup1')
    cmd.save(f"aligned_pymol_{pdb2_filename}.pdb", 'sup2')

    shutil.move(f"aligned_pymol_{pdb1_filename}.pdb", f"{directory}/aligned_pymol_{pdb1_filename}.pdb")
    shutil.move(f"aligned_pymol_{pdb2_filename}.pdb", f"{directory}/aligned_pymol_{pdb2_filename}.pdb")

    # Clean up temporary chain files
    tmp_files = [tmp1, tmp2] + glob.glob("TM*")
    for file in tmp_files:
        if os.path.exists(file):
            os.remove(file)
    
    return f"{directory}/aligned_pymol_{pdb1_filename}.pdb", f"{directory}/aligned_pymol_{pdb2_filename}.pdb"


def read_structure(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_file[:-4], pdb_file)
    return structure

def rename_chain_ids(structure):
    for model in structure:
        for i, chain in enumerate(model):
            chain.id = chr(i + 65)  # ASCII 'A' starts at 65

def write_structure(structure, pdb_file):
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_file)

def seq_similarity(seq1,seq2):
    
    alignments = pairwise2.align.globalxx(seq1, seq2)
    best_alignment = alignments[0]
    aligned_seq1 = best_alignment[0]
    aligned_seq2 = best_alignment[1]
    num_matches = 0
    num_mismatches = 0
    for i in range(len(aligned_seq1)):
        if aligned_seq1[i] == aligned_seq2[i]:
            num_matches += 1
        else:
            num_mismatches += 1
    similarity = num_matches / (num_matches + num_mismatches)
    
    return similarity

def chain_similarity(chain1, chain2):
    if chain1 is None or chain2 is None:
        return 0
    else:
        return seq_similarity(''.join(residue.get_resname() for residue in chain1),
                              ''.join(residue.get_resname() for residue in chain2))

def get_all_chains(structure):
    return [chain for model in structure for chain in model]

def get_chain_order(structure1, structure2):
    chains1 = get_all_chains(structure1)
    chains2 = get_all_chains(structure2)

    num_chains1 = len(chains1)
    num_chains2 = len(chains2)

    # Pad the shorter list with None
    if num_chains1 < num_chains2:
        chains1.extend([None] * (num_chains2 - num_chains1))
    elif num_chains2 < num_chains1:
        chains2.extend([None] * (num_chains1 - num_chains2))

    cost_matrix = -np.array([[chain_similarity(chain1, chain2) for chain1 in chains1]
                             for chain2 in chains2])
    row_ind, col_ind = linear_sum_assignment(cost_matrix)
    return [(chains1[i], chains2[j]) for i, j in zip(row_ind, col_ind)]

def add_chain_to_ligand(pdb_file):
    output_file = f"{pdb_file.split('.')[0]}_ligchain.pdb"
    with open(pdb_file, 'r') as f:
        lines = f.readlines()

    # Find the chain ID
    chain_id = None
    for line in lines:
        if line.startswith("ATOM"):
            chain_id = line[21]  # Chain ID is at position 22 of line (0-indexed)
            break

    # If there is no chain ID, stop the function
    if chain_id is None:
        print("No chain ID found in ATOM records.")
        print(f"Chain ID '{chain_id}' has been added to HETATM records in {output_file}.")
        return pdb_file

    # Add chain ID to HETATM records
    with open(output_file, 'w') as f:
        for line in lines:
            if line.startswith("HETATM") and line[21] == ' ':
                # The chain ID is inserted after the residue name and before the residue number
                newline = line[:21] + chain_id + line[22:]
                f.write(newline)
            else:
                f.write(line)

    return output_file

def detach_hetatms(chain):
    hetatms = []
    residues = list(chain)
    for residue in residues:
        if residue.id[0].startswith("H_"):
            hetatms.append(residue)
            chain.detach_child(residue.id)
    return hetatms

def reattach_hetatms(chain, hetatms):
    for hetatm in hetatms:
        chain.add(hetatm)

def rearrange_TER(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    corrected_lines = []
    atom_flag = False
    hetatm_flag = False
    last_atom_line = None
    for i in range(len(lines)):
        line = lines[i]
        if line.startswith('ATOM'):
            atom_flag = True
            last_atom_line = line
            corrected_lines.append(line)
        elif line.startswith('HETATM'):
            if atom_flag and not hetatm_flag:
                # Generate TER record
                serial_number = int(last_atom_line[6:11].strip()) + 1
                res_name = last_atom_line[17:20]
                chain_id = last_atom_line[21]
                res_seq = last_atom_line[22:26]
                corrected_lines.append(f'TER   {serial_number:5}      {res_name} {chain_id}{res_seq}\n')
                hetatm_flag = True
            corrected_lines.append(line)
        elif line.startswith('TER'):
            # Check if followed by 'ATOM' record
            if i < len(lines) - 1 and lines[i+1].startswith('ATOM'):
                corrected_lines.append(line)
        elif line.startswith('END'):
            corrected_lines.append(line)
        else:
            corrected_lines.append(line)

    with open(filepath, 'w') as f:
        f.writelines(corrected_lines)

def consensus_hets(pdb_file: str, true_pdb: str, big_hets=False):
    pymol.finish_launching(['pymol', '-qc'])  # Quiet and no GUI

    common_ions = ['ZN', 'FE', 'MG', 'MN', 'K', 'NA', 'CA']  # List of ion types to consider

    def get_ions_and_ligand(pdb_file, obj_name):
        cmd.load(pdb_file, obj_name)
        ion_types_present = []

        for ion in common_ions:
            if cmd.count_atoms(f'{obj_name} and resn {ion}'):
                ion_types_present.append(ion)

        ligand_name = ''
        if big_hets:
            ligands = cmd.get_names('objects', 1, f'({obj_name} and hetatm)')
            if ligands:
                largest_ligand = max(ligands, key=lambda name: cmd.count_atoms(f'{name} and hetatm'))
                cmd.set_name(largest_ligand, f'{obj_name}_largest_ligand')
                ligand_name = f'{obj_name}_largest_ligand'

        return ion_types_present, ligand_name

    # Load true_pdb and get ions
    ions2, ligand2 = get_ions_and_ligand(true_pdb, 'molecule2')

    if ions2:
        ion_selections_from_true_pdb = ' or '.join([f'molecule2 and resn {ion}' for ion in ions2])
    else:
        ion_selections_from_true_pdb = ''

    # Process pdb_file
    cmd.load(pdb_file, 'molecule1')
    selections = ['molecule1']
    if ion_selections_from_true_pdb:
        selections.append(ion_selections_from_true_pdb)
    if big_hets:
        ligand1 = get_ions_and_ligand(pdb_file, 'molecule1')[1]
        selections.append(ligand1)
    
    cmd.create('molecule1_with_extras', ' or '.join(selections))
    
    # Save combined molecule1 to a new PDB file
    output1_filename = pdb_file.replace('.pdb', '_modified.pdb')
    cmd.save(output1_filename, 'molecule1_with_extras')

    # Process true_pdb (essentially re-adding the ions and possibly the biggest ligand)
    selections = ['molecule2']
    if ion_selections_from_true_pdb:
        selections.append(ion_selections_from_true_pdb)
    if big_hets:
        selections.append(ligand2)
    
    cmd.create('molecule2_with_extras', ' or '.join(selections))
    
    # Save combined molecule2 to a new PDB file
    output2_filename = true_pdb.replace('.pdb', '_modified.pdb')
    cmd.save(output2_filename, 'molecule2_with_extras')

    # Cleanup
    cmd.delete('molecule1')
    cmd.delete('molecule2')
    cmd.delete('molecule1_with_extras')
    cmd.delete('molecule2_with_extras')

    return output1_filename, output2_filename

def align_and_filter(pdb1, pdb2, big_hets=False):

    pdb1 = add_chain_to_ligand(pdb1)
    pdb2 = add_chain_to_ligand(pdb2)    
   
    pdb1 = uniquify_pdb(pdb1)
    pdb2 = uniquify_pdb(pdb2)
 
    rearrange_TER(pdb1)
    rearrange_TER(pdb2) 

    pdb1, pdb2 = superimpose(pdb1, pdb2)

    pdb1, pdb2 = consensus_hets(pdb1, pdb2, big_hets)

    def residues_to_seq(chain_name):
        """Extracts sequence from the structure in PyMOL and converts to single-letter code."""
        space = {'residues_list': []}
        cmd.iterate(chain_name + ' and name ca', 'residues_list.append(resn)', space=space)
        return ''.join([robust_three_to_one(resn) for resn in space['residues_list']])

    def get_unaligned_residues(seq1, seq2):
        """Returns positions of unaligned residues based on sequence alignment."""
        alignment = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)[0]

        unaligned_positions1 = []
        unaligned_positions2 = []

        resi_seq1 = 0
        resi_seq2 = 0

        for i in range(len(alignment[0])):
            # If there's a gap in seq1's alignment, it means an amino acid from seq2 doesn't align with seq1
            if alignment[0][i] == '-':
                unaligned_positions2.append(resi_seq2 + 1)
                resi_seq2 += 1
            # If there's a gap in seq2's alignment, it means an amino acid from seq1 doesn't align with seq2
            elif alignment[1][i] == '-':
                unaligned_positions1.append(resi_seq1 + 1)
                resi_seq1 += 1
            # If there are no gaps, it means the amino acids from both sequences at the current position align
            else:
                resi_seq1 += 1
                resi_seq2 += 1

        return unaligned_positions1, unaligned_positions2

    def get_chain_sequences(pdb_filename):
        """Get sequences for each chain in a structure using Biopython."""
        chain_sequences = {}
        parser = PDB.PDBParser()
        structure = parser.get_structure(pdb_filename, pdb_filename)

        for model in structure:
            for chain in model:
                seq = []
                for residue in chain:
                    if PDB.is_aa(residue):  # Check if the residue is an amino acid and not a HETATM
                        seq.append(residue.get_resname())
                chain_sequences[chain.id] = ''.join([robust_three_to_one(res) for res in seq])
        return chain_sequences

    def find_best_matching_chain(seq, chain_sequences):
        """Find the best-matching chain based on sequence similarity."""
        best_chain = None
        best_score = -1
        for chain, chain_seq in chain_sequences.items():
            score = pairwise2.align.globalxx(seq, chain_seq, score_only=True)
            if score > best_score:
                best_chain = chain
                best_score = score
        return best_chain

    # Initialize PyMOL
    #pymol.finish_launching()

    # Load the structures
    cmd.load(pdb1, 'structure1')
    cmd.load(pdb2, 'structure2')

    # Extract chain sequences
    struct1_chains = get_chain_sequences(pdb1)
    struct2_chains = get_chain_sequences(pdb2)

    matched_chains_in_structure2 = set()

    for chain1, seq1 in struct1_chains.items():
        if not seq1:
            print(f"Skipping empty sequence for chain {chain1}")
            continue
        best_match_chain = find_best_matching_chain(seq1, struct2_chains)
        if not best_match_chain:
            continue
        matched_chains_in_structure2.add(best_match_chain)
        seq2 = struct2_chains[best_match_chain]

        # Align structure2 chain onto structure1 chain based on sequence alignment
        cmd.align(f'structure2 and chain {best_match_chain}', f'structure1 and chain {chain1}')

        # Get unaligned residue positions for this chain pair
        unaligned_positions1, unaligned_positions2 = get_unaligned_residues(seq1, seq2)

        # Delete non-aligned residues from both structures for this chain pair
        for pos in unaligned_positions1:
            cmd.remove(f'structure1 and chain {chain1} and resi {pos} and not hetatm')

        for pos in unaligned_positions2:
            cmd.remove(f'structure2 and chain {best_match_chain} and resi {pos} and not hetatm')

    # Remove unmatched chains from structure2
    all_chains_in_structure2 = set(struct2_chains.keys())
    unmatched_chains = all_chains_in_structure2 - matched_chains_in_structure2
    for chain in unmatched_chains:
        cmd.remove(f'structure2 and chain {chain} and not hetatm')

    apo = "common_" + pdb1.split('/')[-1]
    holo = "common_" + pdb2.split('/')[-1]
    pdb1 = pdb1.split('/')[0] + '/' + apo
    pdb2 = pdb2.split('/')[0] + '/' + holo
    cmd.save(pdb1, 'structure1')
    cmd.save(pdb2, 'structure2')
    clean_pdb(pdb1, overwrite=True)
    clean_pdb(pdb2, overwrite=True)

    overwrite_pdb_with_pymol(pdb1)
    overwrite_pdb_with_pymol(pdb2)

    #cmd.quit()
    return os.getcwd() + '/' + pdb1, os.getcwd() + '/' + pdb2
                    

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python script.py pdb1 pdb2")
    else:
        pdb1 = sys.argv[1]
        pdb2 = sys.argv[2]
        result = align_and_filter(pdb1, pdb2)
        print("Result:", result)
