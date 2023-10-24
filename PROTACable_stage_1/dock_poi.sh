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
conda activate vs

# Function to calculate centroid of HETATMs
centroid_coordinates() {
    pdb_file="$1"
    awk '
    /^HETATM/ {
        count += 1
        sumx += $7
        sumy += $8
        sumz += $9
    }
    END {
        if (count > 0) {
            print sumx/count, sumy/count, sumz/count
        } else {
            exit 1
        }
    }' "$pdb_file"
}

# Create a temporary PDB without HETATMs
create_protein_pdb() {
    pdb_file="$1"
    tmp_pdb="clean_protein.pdb"
    awk '!/^HETATM/' "$pdb_file" > "$tmp_pdb"
    echo "$tmp_pdb"
}


# Ligand docking function
ligand_docking() {
    echo "Docking the $input set..."
    cp ${PROTACable}/PROTACable_stage_1/utilities/gnina.dpf ./
    cp ${PROTACable}/PROTACable_stage_1/utilities/run_gnina.sh ./

    # modify values in gnina.dpf using sed
    sed -i "s/center_x\s*=\s*/center_x = $CENTER_X/" gnina.dpf
    sed -i "s/center_y\s*=\s*/center_y = $CENTER_Y/" gnina.dpf
    sed -i "s/center_z\s*=\s*/center_z = $CENTER_Z/" gnina.dpf
    sed -i "s/cpu\s*=\s*/cpu = $WORKERS/" gnina.dpf

    out_prefix="docked"

    # modify input and output file names in run_gnina.sh using sed
    sed -i "s/-r\s*protein.pdb/-r ${PROTEIN}/" run_gnina.sh
    sed -i "s/-l\s*library.sdf/-l ${sdf}/" run_gnina.sh
    sed -i "s/-o\s*output.sdf/-o ${out_prefix}_${input%.*}_${sdf}/" run_gnina.sh
    sed -i "s/output.sdf/${out_prefix}_${input%.*}/" run_gnina.sh
    sed -i "s/--cpu\s*cpus/--cpu $WORKERS/" run_gnina.sh
    
    sh run_gnina.sh
    local output="${out_prefix}_${input%.*}_${sdf}"
    obabel -isd $output -opdb -O ${output%.*}.pdb
    echo "$output"
}

# Parsing command-line arguments
while getopts "i:l:x:y:z:ao:w:h" opt; do
    case $opt in
        i) input="$OPTARG";;
        l) sdf="$OPTARG";;
        x) CENTER_X="$OPTARG";;
        y) CENTER_Y="$OPTARG";;
        z) CENTER_Z="$OPTARG";;
        a) AUTO=1;;
        o) OUTPUT="$OPTARG";;
        w) WORKERS="$OPTARG";;
        h)
          echo "Usage: $0 [options]"
          echo "Options:"
          echo "-i  input PDB file. With arbitrary HETATMs to help find the grid if [-a] is provided. Otherwise, without and should provide grid coordinates instead."
          echo "-l  ligand in SDF format."
          echo "-x  grid center x-coordinate."
          echo "-y  grid center y-coordinate."
          echo "-z  grid center z-coordinate."
          echo "-a  automatically infer the grid using the input PDB-ligand complex."
          echo "-o  output docked file."
          echo "-w  number of workers for gnina docking."
          echo "-h  display this help and exit."
          exit 0
          ;;
        *) 
          echo "Invalid option. Use -h for help."
          exit 1
          ;;
    esac
done

if [[ -z "$input" || -z "$sdf" || -z "$WORKERS" ]]; then
    echo "Error: Mandatory arguments -i (input PDB file), -l (ligand in SDF format), and -w (number of workers) must be provided."
    exit 1
fi

# Automatically calculate the centroid if -a is set
if [[ "$AUTO" -eq 1 ]]; then
    read -r CENTER_X CENTER_Y CENTER_Z <<< "$(centroid_coordinates "$input")"
    if [[ $? -ne 0 ]]; then
        echo "Error computing the centroid from the PDB file."
        exit 1
    fi
    PROTEIN=$(create_protein_pdb "$input")
else
    PROTEIN="$input"
fi

echo "Docking ligand with GNINA. Please give it a moment..."
docked=$(ligand_docking)


# If a temporary protein PDB was created, clean it up
if [[ "$AUTO" -eq 1 && -f "$PROTEIN" ]]; then
    rm -f "$PROTEIN"
fi


