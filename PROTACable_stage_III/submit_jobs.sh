#!/bin/bash

#SBATCH --partition=normal
#SBATCH --job-name=P-LINKER

# NEEDS PRESETUP FOR SLURM



source ~/.bashrc
env_exists=$(conda info --envs | grep '^PROTACable ' | wc -l)
if [[ $env_exists -eq 0 ]]; then
    echo "Error: The 'PROTACable' conda environment does not exist."
    exit 1
else
    conda activate PROTACable
fi


pocket=$1
poi_lig=$2
e3_lig=$3


python -W ignore ${PROTACable}/PROTACable_stage_3/main.py $pocket $poi_lig $e3_lig
