"""
------------------------------------------------------------------------------------------
Copyright (C) 2023 Hazem Mslati. All Rights Reserved.

This script is provided "AS IS" without any warranty.
Unauthorized copying, modification, or distribution is prohibited.
------------------------------------------------------------------------------------------
"""

import argparse

def rlabel_pdb(pdb_file_path, atom_index, linker_present):
    with open(pdb_file_path, 'r') as f:
        lines = f.readlines()

    # Find the line with the given atom index
    target_line_index = None
    for i, line in enumerate(lines):
        if line.startswith("HETATM") and int(line[6:11].strip()) == atom_index:
            target_line_index = i
            break

    if target_line_index is None:
        raise ValueError(f"No HETATM line found with atom index {atom_index}.")

    line = lines[target_line_index]

    # Modify atom name based on linker_present value
    modifier = 90 if linker_present else 94
    atom_name = line[12:16].strip()  # Extract the atom name from columns 13 to 16 (0-indexed)
    new_atom_name = f"{atom_name}{modifier}"

    modified_line = f"{line[:12]} {new_atom_name:3s}{line[16:]}"
    
    lines[target_line_index] = modified_line

    new_file_path = pdb_file_path.replace('.pdb', '_radiolabeled.pdb')
    with open(new_file_path, 'w') as f:
        f.writelines(lines)
    
    print(f"Modified pdb file saved to {new_file_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Radiolabel a specified atom in a PDB file based on its atom index.")
    parser.add_argument("pdb_file_path", type=str, help="Path to the input PDB file.")
    parser.add_argument("atom_index", type=int, help="Index of the atom to be modified.")
    parser.add_argument("--linker_present", action='store_true', default=True,
                        help="Specify this flag if a linker is present. Defaults to False if not provided.")
    
    args = parser.parse_args()
    
    rlabel_pdb(args.pdb_file_path, args.atom_index, args.linker_present)
