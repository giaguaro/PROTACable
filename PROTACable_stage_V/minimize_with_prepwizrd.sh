#!/bin/bash

# ------------------------------------------------------------------------------------------
# Copyright (C) 2023 Hazem Mslati. All Rights Reserved.
#
# This script is provided "AS IS" without any warranty.
# Unauthorized copying, modification, or distribution is prohibited.
# ------------------------------------------------------------------------------------------




# Check if the SCHRODINGER environment variable is set
if [ -z "$SCHRODINGER" ]; then
    echo "Error: The SCHRODINGER environment variable is not set."
    exit 1
fi

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <directory>"
    exit 1
fi

DIRECTORY="$1"

if [ ! -d "$DIRECTORY" ]; then
    echo "Error: $DIRECTORY is not a valid directory."
    exit 1
fi

mkdir ternaries_minimzed

for input_protein in "$DIRECTORY"/*.pdb; do
    # Check if the file exists (this avoids issues if there are no .pdb files)
    if [[ -f $input_protein ]]; then
        # Overwrite atom bond templates - for PROTACs
        $SCHRODINGER/utilities/prepwizard "$input_protein" "${input_protein%.*}_pre-minimized.pdb" -NOJOBID -noepik -noccd -f 3 -rmsd 5.0

        $SCHRODINGER/utilities/prepwizard "$input_protein" "${input_protein%.*}_minimized.pdb" -NOJOBID -f 3 -rmsd 5.0
    fi
done

rm ${DIRECTORY}/*pre-minimized.pdb

cp -r ${DIRECTORY}/*_minimized.pdb ternaries_minimized
