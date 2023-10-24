"""
------------------------------------------------------------------------------------------
Copyright (C) 2023 Hazem Mslati. All Rights Reserved.

This script is provided "AS IS" without any warranty.
Unauthorized copying, modification, or distribution is prohibited.
------------------------------------------------------------------------------------------
"""

import pandas as pd
import argparse
import sys
import os
import glob, os
import numpy as np
from biopandas.pdb import PandasPdb


ppdb = PandasPdb()


pdb_df=PandasPdb().read_pdb(str(sys.argv[1]))


pdb_df.df['ATOM']['x_cord']=pdb_df.df['ATOM'].x_coord.round(0)
pdb_df.df['ATOM']['y_cord']=pdb_df.df['ATOM'].y_coord.round(0)
pdb_df.df['ATOM']['z_cord']=pdb_df.df['ATOM'].z_coord.round(0)

df=list(pdb_df.df['ATOM'][pdb_df.df['ATOM'].duplicated(subset=['x_cord','y_cord','z_cord'], keep="first")].atom_number)

pdb_df.df['ATOM'] = pdb_df.df['ATOM'][~pdb_df.df['ATOM'].atom_number.isin(df)]
pdb_df.df['ATOM']=pdb_df.df['ATOM'].drop(columns=['x_cord','y_cord','z_cord'])
pdb_df.to_pdb(f"unclashed_{str(sys.argv[1]).split('.')[0]}.pdb")


