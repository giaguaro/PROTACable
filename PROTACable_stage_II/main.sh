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

UTILS_DIR="${PROTACable}"
H_DEFAULT="${UTILS_DIR}/PROTACable_stage_2/ProPOSE"
L_DEFAULT="${UTILS_DIR}/PROTACable_stage_2/ProPOSE/ProPOSE.lic"

while getopts ":P:H:L:h" opt; do
  case $opt in
    P) P="$OPTARG"
    ;;
    H) H="$OPTARG"
    ;;
    L) L="$OPTARG"
    ;;
    h)
        echo "Usage:"
        echo "-P : Input POI protein-ligand complex to be docked against the E3 Ligase library."
        echo "-H : Where is your ProPOSE software copy installed? Enter HOME directory. Default: ${H_DEFAULT}"
        echo "-L : ProPOSE.lic license path. Default: ${L_DEFAULT}"
        exit 0
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
  esac
done

H="${H:-$H_DEFAULT}"
L="${L:-$L_DEFAULT}"

if [[ -z "$P" ]]; then
    echo "The -P argument (POI protein-ligand complex) must be provided."
    exit 1
fi

id=$(echo "$P" | cut -d '.' -f1 | rev | cut -d '/' -f1 | rev)

for e3 in ${UTILS_DIR}/PROTACable_stage_2/e3_library/*pdb; do
    e=$(basename ${e3}) # same as above ?
    mkdir ${e%.*}_${id}_pair
    cp $e3 ${e%.*}_${id}_pair
    #cp ${UTILS_DIR}/PROTACable_stage_2/ProPOSE* ${e%.*}_${id}_pair
   # cp ${UTILS_DIR}/PROTACable_stage_2/*sh ${e%.*}_${id}_pair
   # cp ${UTILS_DIR}/PROTACable_stage_2/*py ${e%.*}_${id}_pair
    cp $P ${e%.*}_${id}_pair
done

mkdir ternaries

for pair in *_pair; do
    (cd $pair && sbatch ${UTILS_DIR}/PROTACable_stage_2/bash_propose.sh -p $P -e e3*pdb -H $H -L $L)
done

