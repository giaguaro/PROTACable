"""
------------------------------------------------------------------------------------------
Copyright (C) 2023 Hazem Mslati. All Rights Reserved.

This script is provided "AS IS" without any warranty.
Unauthorized copying, modification, or distribution is prohibited.
------------------------------------------------------------------------------------------
"""

import pandas as pd
import glob, os
from biopandas.mol2 import PandasMol2
from biopandas.mol2 import split_multimol2
import numpy as np
from scipy.spatial import distance

pdmol = PandasMol2()

def save_e3_mols (file, e3_glob):
    with open(f'{file}_filtered.mol2', 'w') as f:
        for emol2 in e3_glob:
            e3_pd=PandasMol2().read_mol2_from_list(mol2_lines=emol2[1], mol2_code=emol2[0])
            f.write(e3_pd.mol2_text)
        
        
#for pp in glob.glob('./*_pair'):
#    print(pp)
e3_glob=glob.glob(f'./e3_*p1fixed.mol2')[0]
e3=glob.glob(f'./e3_*p1fixed.mol2')[0].split('/')[-1].split('.')[0]
poi_glob=glob.glob(f'./poi_*p1fixed.mol2')[0]
poi=glob.glob(f'./poi_*p1fixed.mol2')[0].split('/')[-1].split('.')[0]
df_scores=pd.read_csv(f'./pp_scores.csv', sep=',')
try:
    os.remove(f'./{poi}_filtered.mol2')
except:
    pass
try: 
    os.remove(f'./{e3}_filtered.mol2')
except:
    pass
try:
    os.remove(f'./pp_scores_filtered.csv')
except:
    pass
with open(f'./{poi}_filtered.mol2', 'w') as f:
    for pmol2 in split_multimol2(f'{poi_glob}'):
        poi_pd=PandasMol2().read_mol2_from_list(mol2_lines=pmol2[1], mol2_code=pmol2[0])
        f.write(poi_pd.mol2_text)
        pligand_coor=poi_pd.df[poi_pd.df['subst_name'].str.contains("UNL")].iloc[:,2:5].values
        with open(f'./{e3}_filtered.mol2', 'w') as f:
            print(poi)
            list_e3=[]
            #keys=['cityblock']
            #distanc_lis={key: [] for key in keys}
           #print(distanc_lis)
            distanc_lis=[]

            for emol2 in split_multimol2(f'{e3_glob}'):
                    keep_molecule=False
                    e3_pd=PandasMol2().read_mol2_from_list(mol2_lines=emol2[1], mol2_code=emol2[0])

                    # do some analysis

                    keto_coord = e3_pd.df[e3_pd.df['atom_name'] == 'O95'][['x', 'y', 'z']]
                    if keto_coord.empty:
                        keto_coord = e3_pd.df[e3_pd.df['atom_name'] == 'N95'][['x', 'y', 'z']]
                    if keto_coord.empty:
                        keto_coord = e3_pd.df[e3_pd.df['atom_name'] == 'C95'][['x', 'y', 'z']]
                    if keto_coord.empty:
                        keto_coord = e3_pd.df[e3_pd.df['atom_name'] == 'S95'][['x', 'y', 'z']]
                    if keto_coord.empty:
                        keto_coord = e3_pd.df[e3_pd.df['atom_name'] == 'P95'][['x', 'y', 'z']]

                    distances = poi_pd.distance(keto_coord.values[0])
                    poi_pd.df['distances'] = distances

                    dis99=poi_pd.df[poi_pd.df['atom_name'].str.contains("O94|N94|C94|S94|P94")].distances.values


                    eligand_coor=e3_pd.df[e3_pd.df['subst_name'].str.contains("UNL")].iloc[:,2:5].values
                    #possibilities=["braycurtis", "canberra", "chebyshev", 'cityblock', "cosine", "euclidean", "mahalanobis", "minkowski", "seuclidean", 'sqeuclidean']
                    distanc_lis.append((dis99+np.mean(distance.cdist(pligand_coor, eligand_coor, "cityblock"))/2))
                    #distanc_lis.append(dis99)
                        ##distance_centroid=distance.cdist(pligand_coor, eligand_coor, 'euclidean')
                    #print(pligand_coor.shape, eligand_coor.shape,distance_centroid.shape)
                    #distance_centroid=np.linalg.norm(pligand_coor - eligand_coor)
                    #print(distance_centroid)
                    ##distanc_lis.append(min(min(distance_centroid, key=min)))
                    #print(min(min(distance_centroid, key=min)))


#                         if 0<float(dis.item())<float(i):
#                             #print(dis.item())
#                             keep_molecule=True
#                             list_e3.append(e3)
#                             #print(e3)








                    # save molecule if it passes our filter criterion
                    #if keep_molecule: 
                        # note that the mol2_text contains the original mol2 content

                        #print(' do nothing')

                        #f.write(e3_pd.mol2_text) 
#             if not list_e3:
#                 print(e3, " does not have a distance")
#                 print(min(distanc_lis))
        #print(distanc_lis)
        try:
            a=np.asarray(distanc_lis).reshape(1,100)
        except:
            a=np.asarray(distanc_lis)[:,1:].reshape(1,100)
        ind=np.argpartition(a[0], 19)[:20].tolist()
        print (ind)
        multimol2=[list(split_multimol2(f'{e3_glob}'))[i] for i in ind]
        file=f'./{e3}'
        save_e3_mols(file,multimol2)
        pp_scores=df_scores.iloc[ind]
        pp_scores.to_csv(f'./pp_scores_filtered.csv', index=False)
        #print(pp)

        break
