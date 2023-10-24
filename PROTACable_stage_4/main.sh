#!/bin/bash

# ------------------------------------------------------------------------------------------
# Copyright (C) 2023 Hazem Mslati. All Rights Reserved.
#
# This script is provided "AS IS" without any warranty.
# Unauthorized copying, modification, or distribution is prohibited.
# ------------------------------------------------------------------------------------------


source ~/.bashrc
conda activate PROTACable

if [[ -z "${PROTACable}" ]]; then
    echo "Please export the PROTACable variable to point to the PROTACable home script directory."
    exit 1
fi

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_directory>"
    exit 1
fi

input_dir_="$1"
input_dir=$(realpath "$input_dir_")
ligand_dir_=${input_dir}/ligands
ligand_dir=$(realpath "$ligand_dir_")
utils_dir="${PROTACable}/PROTACable_stage_4/model/embed" 
checkpoint="${PROTACable}/PROTACable_stage_4/model/checkpoint/epoch_50-val_auc_roc_0.936-val_acc_0.913-val_f1_0.923-val_loss_0.302.ckpt"

init_pwd=$(pwd)

for pdb_file in "${input_dir}"/*.pdb; do
    base_name=$(basename "${pdb_file%.*}")
    echo $base_name
    temp_dir="${input_dir}/${base_name}_temp_embedding_dir"

    mkdir -p "${temp_dir}"

    cp "${pdb_file}" "${temp_dir}/"
    cd "${temp_dir}/"

    python "${utils_dir}/process_datapoint.py" "${base_name}.pdb" --allow_exceptions True --protacs_embedding True --ligand_docked True --clean True --interface 400 --featurize True --inference True --ligands_path $ligand_dir --by_id $base_name

    pkl_file=$(ls *.pkl)
    if [ -f "${pkl_file}" ]; then
        cp "${pkl_file}" "${input_dir}/"
    fi

    #rm -r "${temp_dir}"
    cd $init_pwd
done

for pkl_file in "${input_dir}"/*.pkl; do
    echo "${pkl_file}" >> "${input_dir}/raw_prediction_paths.txt"
done

echo "Embedding computation completed!"


echo "Predicting ..."

cd ${input_dir}
python -u ${utils_dir}/protac_predict.py --checkpoint_path "${checkpoint}" --test_file "${input_dir}/raw_prediction_paths.txt"


rm -r *pkl *dir


mkdir ../top_ternaries_results/

awk -F',' 'NR>1 {print $2, $1}' ${input_dir}/predictions_and_targets.csv | sort -nr | head -20 | awk '{print $2}' | sed 's/_clean_interface_lframe.pdb/.pdb/' | xargs -I{} cp {} ../top_ternaries_results/

