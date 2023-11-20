#!/bin/bash

#SBATCH --job-name=GNINA
#SBATCH --mem=0

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./extra_libs/
${PROTACable}/PROTACable_stage_1/utilities/gnina --cpu cpus --no_gpu --cnn_scoring none -r protein.pdb -l library.sdf -o output.sdf --config gnina.dpf


