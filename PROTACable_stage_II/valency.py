"""
------------------------------------------------------------------------------------------
Copyright (C) 2023 Hazem Mslati. All Rights Reserved.

This script is provided "AS IS" without any warranty.
Unauthorized copying, modification, or distribution is prohibited.
------------------------------------------------------------------------------------------
"""

import pandas as pd
import glob, os
import sys
from biopandas.pdb import PandasPdb


pdb_df=PandasPdb().read_pdb(str(sys.argv[1]))

#pdb_df=PandasPdb().read_pdb('/groups/cherkasvgrp/share/progressive_docking/hmslati/protacs/manual_linker_disjoint_propose_docking/maes/mae_amb_complex_poi_7nq8_9_docked_e3_cereblon_7nq8_9_docked_pp_model_xx14_mae.pdb')
reference_pointP = pdb_df.df['HETATM'][pdb_df.df['HETATM']['atom_name'].str.contains("O94|N94|C94|S94|P94")][['x_coord', 'y_coord', 'z_coord']].to_numpy()[0]
distances = pdb_df.distance(xyz=reference_pointP, records=('HETATM',))
all_within_1pt5A0 = pdb_df.df['HETATM'][distances < 1.5]
deli0=all_within_1pt5A0[all_within_1pt5A0['atom_name'].str.contains("H")].iloc[0].atom_number


reference_pointE = pdb_df.df['HETATM'][pdb_df.df['HETATM']['atom_name'].str.contains("O95|N95|C95|S95|P95")][['x_coord', 'y_coord', 'z_coord']].to_numpy()[0]
distances = pdb_df.distance(xyz=reference_pointE, records=('HETATM',))
all_within_1pt5A1 = pdb_df.df['HETATM'][distances < 1.5]
deli1=all_within_1pt5A1[all_within_1pt5A1['atom_name'].str.contains("H")].iloc[0].atom_number
delis=[deli0,deli1]
pdb_df.df['HETATM'] = pdb_df.df['HETATM'][~pdb_df.df['HETATM'].atom_number.isin(delis)]

pdb_df.to_pdb(f"valency_{str(sys.argv[1]).split('.')[0]}.pdb")
