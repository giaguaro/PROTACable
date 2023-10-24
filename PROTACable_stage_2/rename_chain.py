"""
------------------------------------------------------------------------------------------
Copyright (C) 2023 Hazem Mslati. All Rights Reserved.

This script is provided "AS IS" without any warranty.
Unauthorized copying, modification, or distribution is prohibited.
------------------------------------------------------------------------------------------
"""

import sys
from biopandas.pdb import PandasPdb
import argparse

def rename_chain(pdb_file, chain_name):
    ppdb = PandasPdb().read_pdb(pdb_file)

    for df_name in ppdb.df.keys():
        if 'chain_id' in ppdb.df[df_name].columns:
            ppdb.df[df_name]['chain_id'] = chain_name

    ppdb.to_pdb(path=pdb_file, gz=False, append_newline=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Rename chain in a PDB file.')
    parser.add_argument('pdb_file', type=str, help='Input PDB file')
    parser.add_argument('chain_name', type=str, help='New chain name')
    
    args = parser.parse_args()

    if len(args.chain_name) != 1:
        print("Error: Chain name should be a single character.")
        sys.exit(1)

    rename_chain(args.pdb_file, args.chain_name)
