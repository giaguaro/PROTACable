#!/bin/bash

# ------------------------------------------------------------------------------------------
# Copyright (C) 2023 Hazem Mslati. All Rights Reserved.
#
# This script is provided "AS IS" without any warranty.
# Unauthorized copying, modification, or distribution is prohibited.
# ------------------------------------------------------------------------------------------


##SBATCH --partition=<
#SBATCH --job-name=LEAP

source ~/.bashrc
conda activate PROTACable

p=$1
e=$2

cp ${p%.*}_r_prep_pp_p1fixed_filtered.mol2 ${p%.*}_r_prep_pp_fixed.mol2
cp ${e%.*}_l_prep_pp_p1fixed_filtered.mol2 ${e%.*}_l_prep_pp_fixed.mol2

poi="${p%.*}_r_prep_pp_fixed.mol2"
e3="${e%.*}_l_prep_pp_fixed.mol2"

rm xx* *model* diagnosis* *temp*
csplit $e3 '/@<TRIPOS>MOLECULE/' '{*}'
rm xx00

obabel -imol2 $poi -opdb -O ${poi%.*}.pdb
sed -i -e '/UNL/ s/ATOM  /HETATM/' ${poi%.*}.pdb
sed -i '/HETATM/d' ${poi%.*}.pdb
sed -i '/CONECT/d' ${poi%.*}.pdb
sed -i '/MASTER/d' ${poi%.*}.pdb
sed -i '/END/d' ${poi%.*}.pdb
${PROTACable}/PROTACable_stage_2/rename_chain.py ${poi%.*}.pdb "A"

v30=$(echo "$e3" | cut -d '.' -f1 | rev | cut -d '_' -f7-8 | rev)
v310=$(basename $poi)
v31=${v310%.*}
v3=${v30}_${v31}

mkdir ../ternaries/${v3}
mkdir ../ternaries/${v3}/top_20_pooled_${v3}
mkdir ../ternaries/${v3}/top_20_pooled_${v3}/for_linker

grep "^HETATM" $p > ${poi%.*}_temp_lig.pdb
${PROTACable}/PROTACable_stage_2/rename_chain.py ${poi%.*}_temp_lig.pdb "A"

for x in xx*; do
    obabel -imol2 $x -opdb -O ${x%.*}.pdb
    sed -i -e '/UNL/ s/ATOM  /HETATM/' ${x%.*}.pdb
    sed -i '/HETATM/d' ${x%.*}.pdb
    ${PROTACable}/PROTACable_stage_2/align.py $e ${x%.*}.pdb
    ${PROTACable}/PROTACable_stage_2/rename_chain.py aligned_complex.pdb "B"
    grep "^HETATM" aligned_complex.pdb >> temp_lig2.pdb
    sed -i '/HETATM\|CONECT\|END/d' aligned_complex.pdb

    cat ${poi%.*}.pdb ${poi%.*}_temp_lig.pdb >> complex_${v3}_pp_model_${x%.*}.pdb
    echo "TER" >>  complex_${v3}_pp_model_${x%.*}.pdb
    cat aligned_complex.pdb temp_lig2.pdb >> complex_${v3}_pp_model_${x%.*}.pdb
    echo "END" >> complex_${v3}_pp_model_${x%.*}.pdb

    pdb4amber -i complex_${v3}_pp_model_${x%.*}.pdb -o amb_complex_${v3}_pp_model_${x%.*}.pdb
    mv amb_complex_${v3}_pp_model_${x%.*}.pdb complex_${v3}_pp_model_${x%.*}.pdb
    cp complex_${v3}_pp_model_${x%.*}.pdb ../ternaries/${v3}/top_20_pooled_${v3}

    grep "^ATOM" complex_${v3}_pp_model_${x%.*}.pdb >> pocket_pp_model_${x%.*}.pdb
    cp pocket_pp_model_${x%.*}.pdb ../ternaries/${v3}/top_20_pooled_${v3}/for_linker
    cp ${poi%.*}_temp_lig.pdb ../ternaries/${v3}/top_20_pooled_${v3}/for_linker/${poi%.*}_poi_lig_pp_model.pdb
    cp temp_lig2.pdb ../ternaries/${v3}/top_20_pooled_${v3}/for_linker/e3_lig_pp_model_${x%.*}.pdb

    rm -r ${x%.*}.pdb temp_lig2.pdb *log *XX*
done

rm ${poi%.*}_temp_lig.pdb *xx* *temp*

echo "Done. See results in ternaries/${v3}/top_20_pooled_${v3} directory"
