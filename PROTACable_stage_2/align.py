"""
------------------------------------------------------------------------------------------
Copyright (C) 2023 Hazem Mslati. All Rights Reserved.

This script is provided "AS IS" without any warranty.
Unauthorized copying, modification, or distribution is prohibited.
------------------------------------------------------------------------------------------
"""

import sys
from pymol import cmd, stored

def align_complex_to_protein(mobile_complex, immobile_protein):
    """
    Aligns the protein part of a protein-ligand complex to a reference protein.
    The ligand's pose in the protein-ligand complex will be preserved.
    """
    cmd.load(mobile_complex, 'mobile')
    cmd.load(immobile_protein, 'immobile')
    
    cmd.select('mobile_protein', 'mobile and polymer')
    cmd.select('mobile_ligand', 'mobile and not polymer')
    
    cmd.align('mobile_protein', 'immobile')
    
    cmd.save('aligned_complex.pdb', 'mobile')

    print("Alignment complete. The aligned complex has been saved as 'aligned_complex.pdb'.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python align_script.py <path_to_mobile_complex> <path_to_immobile_protein>")
        sys.exit(1)

    mobile_complex_path = sys.argv[1]
    immobile_protein_path = sys.argv[2]
    
    align_complex_to_protein(mobile_complex_path, immobile_protein_path)

