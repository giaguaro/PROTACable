"""
------------------------------------------------------------------------------------------
Copyright (C) 2023 Hazem Mslati. All Rights Reserved.

This script is provided "AS IS" without any warranty.
Unauthorized copying, modification, or distribution is prohibited.
------------------------------------------------------------------------------------------
"""

##selects the shortest distance between 90 (POI) & 91 (Linker) attachment points##

##############################################################################################################
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
import glob
import numpy as np
from biopandas.pdb import PandasPdb


for sd in glob.glob(f'./linker_library/temp/*docked*.sdf'):
    name=sd.split('/')[-1].split('.')[0]
    cname='_'.join(name.split('_')[0:2])
    print(cname)
    for c in glob.glob(f'./conect_cp/{cname}_cnct.idx'):
        distances={}
        infile = open(c, 'r')
        firstLine = infile.readline()
        position=firstLine.split(',')[1]
        position0=firstLine.split(',')[0]
        if position:
            mols = Chem.SDMolSupplier(sd)
            for pdb in glob.glob(f'./*_rec.pdb'):
                poi=PandasPdb().read_pdb(f'{pdb}')
                keto_coord = poi.df['HETATM'].iloc[[position0]][['x_coord', 'y_coord', 'z_coord']].to_numpy()
                w = Chem.SDWriter(f'./linker_library/temp/{name}_short.sdf')
                for idx, m in enumerate(mols):
                    pos = m.GetConformer().GetAtomPosition(int(position))
                    one=pos.x, pos.y, pos.z
                    distance = np.sqrt(np.sum(np.square(keto_coord - one)))
                    distances[idx]=distance
                w.write(mols[min(distances, key=distances.get)])
        else:
            pass
        
##############################################################################################################


# DUMMY LINKERS

for sd in glob.glob(f'./dummy_linkers/temp/*docked*.sdf'):
    name=sd.split('/')[-1].split('.')[0]
    cname='_'.join(name.split('_')[0:2])
    print(cname)
    for c in glob.glob(f'./dummy_conect/{cname}_cnct.idx'):
        distances={}
        infile = open(c, 'r')
        firstLine = infile.readline()
        position=firstLine.split(',')[1]
        position0=firstLine.split(',')[0]
        if position:
            mols = Chem.SDMolSupplier(sd)
            for pdb in glob.glob(f'./*_rec.pdb'):
                poi=PandasPdb().read_pdb(f'{pdb}')
                keto_coord = poi.df['HETATM'].iloc[[position0]][['x_coord', 'y_coord', 'z_coord']].to_numpy()
                w = Chem.SDWriter(f'./dummy_linkers/temp/{name}_short.sdf')
                for idx, m in enumerate(mols):
                    pos = m.GetConformer().GetAtomPosition(int(position))
                    one=pos.x, pos.y, pos.z
                    distance = np.sqrt(np.sum(np.square(keto_coord - one)))
                    distances[idx]=distance
                w.write(mols[min(distances, key=distances.get)])
        else:
            pass
