"""
------------------------------------------------------------------------------------------
Copyright (C) 2023 Hazem Mslati. All Rights Reserved.

This script is provided "AS IS" without any warranty.
Unauthorized copying, modification, or distribution is prohibited.
------------------------------------------------------------------------------------------
"""

import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem

def show_3d_mol(mol, width=300, height=300):
#     mol = Chem.AddHs(mol)  # Add hydrogens to fill valence
#     AllChem.EmbedMolecule(mol, AllChem.ETKDG())  
    block = Chem.MolToMolBlock(mol)
    viewer = py3Dmol.view(width=width, height=height)
    viewer.addModel(block, format="mol")
    viewer.setStyle({'stick': {}})
    viewer.zoomTo()
    return viewer.show()
