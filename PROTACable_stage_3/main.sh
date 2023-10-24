#!/bin/bash

# ------------------------------------------------------------------------------------------
# Copyright (C) 2023 Hazem Mslati. All Rights Reserved.
#
# This script is provided "AS IS" without any warranty.
# Unauthorized copying, modification, or distribution is prohibited.
# ------------------------------------------------------------------------------------------

if [[ -z "${PROTACable}" ]]; then
    echo "Please export the PROTACable variable to point to the PROTACable home script directory."
    exit 1
fi


source ~/.bashrc
env_exists=$(conda info --envs | grep '^PROTACable ' | wc -l)
if [[ $env_exists -eq 0 ]]; then
    echo "Error: The 'amber' conda environment does not exist."
    exit 1
else
    conda activate PROTACable
fi

init_dir=$1

desired_dir="${init_dir%ternaries*}"

cd "$desired_dir" || exit 1

for filepath in ternaries/*/top_20_pooled*/for_linker; do
    poi_lig=$(ls ${filepath}/*poi_lig_pp_model.pdb)
    for pocket in ${filepath}/pocket_pp_model*pdb; do
        pocket_bname=$(basename ${pocket})
        model_=$(echo ${pocket_bname}| cut -d'_' -f4)
        model_nb=${model_%%.*}
        echo "Now processing ${filepath} with model number ${model_nb}"
        e3_lig=$(ls ${filepath}/e3_lig_pp_model_*${model_nb}.pdb)
        python -W ignore ${PROTACable}/PROTACable_stage_3/main.py $pocket $poi_lig $e3_lig
    done
done

