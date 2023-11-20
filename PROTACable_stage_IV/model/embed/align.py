"""
------------------------------------------------------------------------------------------
Copyright (C) 2023 Hazem Mslati. All Rights Reserved.

This script is provided "AS IS" without any warranty.
Unauthorized copying, modification, or distribution is prohibited.
------------------------------------------------------------------------------------------
"""

import sys
from pymol import cmd

def superimpose_proteins(pdb1, pdb2):
    # Load the two PDBs into PyMOL
    cmd.load(pdb1, "protein1")
    cmd.load(pdb2, "protein2")

    # Remove non-protein atoms
    cmd.remove("protein1 and not polymer.protein")
    cmd.remove("protein2 and not polymer.protein")

    # Superimpose protein2 onto protein1
    cmd.align("protein2 and name CA", "protein1 and name CA")

    # Save the superimposed structures
    cmd.save(pdb1.replace(".pdb", "_super.pdb"), "protein1")
    cmd.save(pdb2.replace(".pdb", "_super.pdb"), "protein2")

    # Close the loaded PDBs
    cmd.delete("protein1")
    cmd.delete("protein2")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python superimpose.py <pdb1_path> <pdb2_path>")
        sys.exit(1)

    pdb1_path = sys.argv[1]
    pdb2_path = sys.argv[2]

    superimpose_proteins(pdb1_path, pdb2_path)

