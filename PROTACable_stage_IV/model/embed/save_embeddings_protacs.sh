#/bin/bash

source ~/.bashrc
conda activate gmx


utils_dir=""

mydir=$1
#cd protacs/ternaries/${mydir}

input=$(ls *-2_complex.pdb)
#output=${input%.*}_interface.pdb

#part=$(echo $input | cut -d'_' -f6-7)
#rm *pkl
python ${utils_dir}process_datapoint.py "$input" --label True --label_file ${utils_dir}protac_activity_db.csv --label_column flag --by_id $part --allow_exceptions True --id_column identifier --protacs_embedding True --ligand_docked True --clean True --interface 400

