"""
------------------------------------------------------------------------------------------
Copyright (C) 2023 Hazem Mslati. All Rights Reserved.

This script is provided "AS IS" without any warranty.
Unauthorized copying, modification, or distribution is prohibited.
------------------------------------------------------------------------------------------
"""

import random
import string
import sys

def random_string(length):
    return ''.join(random.choices(string.ascii_uppercase, k=length))

def rename_residues(pdb_filename):
    new_filename = f"{pdb_filename.split('.')[0]}_renamed.pdb"
    name = random_string(3)
    with open(pdb_filename, 'r') as file:
        lines = file.readlines()
    with open(new_filename, 'w') as file:
        for line in lines:
            if line.startswith("HETATM"):
                line_list = list(line)
                line_list[17:20] = name
                file.write(''.join(line_list))
            else:
                file.write(line)
    return new_filename

if __name__ == '__main__':
    pdb_filename = sys.argv[1]
    rename_residues(pdb_filename)

