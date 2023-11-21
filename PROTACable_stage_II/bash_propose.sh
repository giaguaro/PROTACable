#!/bin/bash

# ------------------------------------------------------------------------------------------
# Copyright (C) 2023 Hazem Mslati. All Rights Reserved.
#
# This script is provided "AS IS" without any warranty.
# Unauthorized copying, modification, or distribution is prohibited.
# ------------------------------------------------------------------------------------------

#SBATCH --partition=normal
#SBATCH --job-name=P-PP_dock


if [[ -z $SLURM_JOB_ID ]]; then
    echo "Error: This script should be run with SLURM."
    exit 1
fi

##SBATCH --partition=<REPLACE_WITH_YOUR_PARTIION>
#SBATCH --job-name=PP_dock

source ~/.bashrc
env_exists=$(conda info --envs | grep '^PROTACable ' | wc -l)
if [[ $env_exists -eq 0 ]]; then
    echo "Error: The 'amber' conda environment does not exist."
    exit 1
else
    conda activate PROTACable
fi

if [[ -z "${PROTACable}" ]]; then
    echo "Please export the PROTACable variable to point to the PROTACable home script directory."
    exit 1
fi

PROPOSEHOME="$ProPOSE"
LICENSE_PATH="${PROPOSEHOME}/ProPOSE.lic"

while getopts ":p:e:H:L:" opt; do
    case $opt in
        p)
            p="$OPTARG"
            ;;
        e)
            e3="$OPTARG"
            ;;
        H)
            PROPOSEHOME="$OPTARG"
            ;;
        L)
            LICENSE_PATH="$OPTARG"
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done

if [[ -z $p || -z $e3 ]]; then
    echo "Error: Both -p and -e arguments are mandatory."
    exit 1
fi

echo "%{j}.out"

# Extract ID from file name
#id=$(echo "$p" | cut -d '.' -f1 | rev | cut -d '_' -f4-5 | rev)
id=$(echo "$p" | cut -d '.' -f1 | rev | cut -d '/' -f1 | rev)

cp "$p" "poi_target_${id}.pdb"
poi="poi_target_${id}.pdb"

# Verify that LICENSE_PATH is provided
if [[ -z $LICENSE_PATH ]]; then
    echo "Error: License path (-L) is mandatory."
    exit 1
fi


export PROPOSEHOME

rm -f *scores.csv *special* *hit *trial* *prep* *ligand* sqm.* ANTECHAMBER* ATOMTYPE* *log *noHET* *sed* *mol2 *sslink* *txt

cp $poi ${poi%.*}_prep1.pdb
cp $e3 ${e3%.*}_prep1.pdb

echo "renaming HETATMs"

rename_atoms() {
    local file_prefix=$1
    local output_file="${file_prefix}_prep1.pdb"
    
    for i in {1..80}; do
        for atm in $o_atoms; do
            local spacing=" "
            [ $i -lt 10 ] && spacing="  "
            sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i}${spacing}UNL/}" $output_file
        done
        for atm in $t_atoms; do
            local spacing=" "
            [ $i -lt 10 ] && spacing="  "
            sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i}${spacing}UNL/}" $output_file
        done
        local spacing=" "
        [ $i -lt 10 ] && spacing="  "
        sed -i "$((i-1)),/ H   UNL/{s/ H   UNL/ H${i}${spacing}UNL/}" $output_file
    done

    sed -i -e 's/ H  1 UNL/  H   UNL/g' -e 's/ H  2 UNL/  H   UNL/g' -e 's/ H  3 UNL/  H   UNL/g' $output_file

    for hyd in $t_Hatoms $th_Hatoms $f_Hatoms; do
        local replacement=" H   UNL"
        [ "${#hyd}" -eq 4 ] && replacement="  H   UNL"
        sed -i "s/ $hyd UNL/$replacement/g" $output_file
    done
}

o_atoms="C O N F S P"
t_atoms="Cl Br"
t_Hatoms="HN HO HS HP"
th_Hatoms="HXT HN1 HN2 HN3 HO1 1HO HO2 2HO HO3 3HO HS1 1HS HS2 2SH HS3 3SH HP1 1HP HP2 2HP HP3 3HP"
f_Hatoms="1HXT 2HXT 3HXT"

rename_atoms "${poi%.*}"
rename_atoms "${e3%.*}"


####################################################################################################################################################

#### spesh case
	
process_special_case() {
    local base=$1
    
    grep "^HETATM" ${base}_prep1_special.pdb > ${base}_prep1_special_lig.pdb
    grep "^CONECT" ${base}_prep1_special.pdb >> ${base}_prep1_special_lig.pdb
    grep "^END" ${base}_prep1_special.pdb >> ${base}_prep1_special_lig.pdb

    for pattern in HETATM CONECT END; do
        sed -i "/$pattern/d" ${base}_prep1_special.pdb
    done

    pdb4amber -i ${base}_prep1_special.pdb -o ${base}_r_special.pdb --nohyd
    
    for pattern in HETATM CONECT END; do
        sed -i "/$pattern/d" ${base}_r_special.pdb
    done

    cat ${base}_prep1_special_lig.pdb >> ${base}_r_special.pdb
    pdb4amber -i ${base}_r_special.pdb -o ${base}_r_special2.pdb
}

# Special Case
cp ${poi%.*}_prep1.pdb ${poi%.*}_prep1_special.pdb
process_special_case ${poi%.*}

# Main Operation
pdb4amber -i ${poi%.*}_prep1.pdb -o ${poi%.*}_r_prep1.pdb

if [ ! -f ${poi%.*}_r_prep1.pdb ] || [ ! -s ${poi%.*}_r_prep1.pdb ]; then
    conda deactivate amber
    conda activate py39

    for pattern in HETATM CONECT END; do
        grep "^$pattern" ${poi%.*}_prep1.pdb >> ${poi%.*}_ligand0.pdb
        sed -i "/$pattern/d" ${poi%.*}_prep1.pdb
    done

    pdbfixer ${poi%.*}_prep1.pdb --replace-nonstandard --add-residues --output="fixed_${poi%.*}_prep1.pdb"
    sed -i '/END/d' fixed_${poi%.*}_prep1.pdb
    sed -i '/TER.*UNL\|UNL.*TER/d' fixed_${poi%.*}_prep1.pdb
    cat ${poi%.*}_ligand0.pdb >> fixed_${poi%.*}_prep1.pdb

    conda deactivate py39
    conda activate amber

    pdb4amber -i fixed_${poi%.*}_prep1.pdb -o ${poi%.*}_r_prep1.pdb
    cp $poi ${poi%.*}_prep1.pdb

    echo "renaming HETATMs"

    o_atoms="C O N F S P"
    t_atoms="Cl Br"
    o_Hatoms="H"
    t_Hatoms="HN HO HS HP"
    th_Hatoms="HXT HN1 HN2 HN3 HO1 1HO HO2 2HO HO3 3HO HS1 1HS HS2 2SH HS3 3SH HP1 1HP HP2 2HP HP3 3HP"
    f_Hatoms="1HXT 2HXT 3HXT"

    #POI
    ####
    for i in {1..9}; do for atm in $o_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i}  UNL/}" ${poi%.*}_prep1.pdb;done;done
    for i in {10..80}; do for atm in $o_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i} UNL/}" ${poi%.*}_prep1.pdb;done;done

    for i in {1..9}; do for atm in $t_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i}  UNL/}" ${poi%.*}_prep1.pdb;done;done
    for i in {10..80}; do for atm in $t_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i} UNL/}" ${poi%.*}_prep1.pdb;done;done

    sed -i 's/ H  1 UNL/  H   UNL/g' ${poi%.*}_prep1.pdb; sed -i 's/ H  2 UNL/  H   UNL/g' ${poi%.*}_prep1.pdb; sed -i 's/ H  3 UNL/  H   UNL/g' ${poi%.*}_prep1.pdb;

    for hyd in $t_Hatoms; do sed -i "s/ $hyd  UNL/ H   UNL/g" ${poi%.*}_prep1.pdb;done
    for hyd in $th_Hatoms; do sed -i "s/ $hyd UNL/ H   UNL/g" ${poi%.*}_prep1.pdb;done
    for hyd in $f_Hatoms; do sed -i "s/ $hyd UNL/  H   UNL/g" ${poi%.*}_prep1.pdb;done


    for i in {1..9}; do sed -i "$((i-1)),/ H   UNL/{s/ H   UNL/ H${i}  UNL/}" ${poi%.*}_prep1.pdb;done
    for i in {10..80}; do sed -i "$((i-1)),/ H   UNL/{s/ H   UNL/ H${i} UNL/}" ${poi%.*}_prep1.pdb;done

    cp fixed_${poi%.*}_prep1.pdb ${poi%.*}_prep1_special.pdb
    rm ${poi%.*}_r_special.pdb ${poi%.*}_r_special2.pdb
    process_special_case ${poi%.*}
fi


pdb4amber -i ${e3%.*}_prep1.pdb -o ${e3%.*}_r_prep1.pdb



#### spesh case
cp ${e3%.*}_prep1.pdb ${e3%.*}_prep1_special.pdb
grep "^HETATM" ${e3%.*}_prep1_special.pdb > ${e3%.*}_prep1_special_lig.pdb
grep "^CONECT" ${e3%.*}_prep1_special.pdb >> ${e3%.*}_prep1_special_lig.pdb
grep "^END" ${e3%.*}_prep1_special.pdb >> ${e3%.*}_prep1_special_lig.pdb

sed -i '/HETATM/d' ${e3%.*}_prep1_special.pdb
sed -i '/CONECT/d' ${e3%.*}_prep1_special.pdb
sed -i '/END/d' ${e3%.*}_prep1_special.pdb


pdb4amber -i ${e3%.*}_prep1_special.pdb -o ${e3%.*}_r_special.pdb --nohyd
sed -i '/HETATM/d' ${e3%.*}_r_special.pdb
sed -i '/CONECT/d' ${e3%.*}_r_special.pdb
sed -i '/END/d' ${e3%.*}_r_special.pdb

cat ${e3%.*}_prep1_special_lig.pdb >> ${e3%.*}_r_special.pdb

pdb4amber -i ${e3%.*}_r_special.pdb -o ${e3%.*}_r_special2.pdb
#####

if [ ! -f ${e3%.*}_r_prep1.pdb ] || [ ! -s ${e3%.*}_r_prep1.pdb ]; then
        conda deactivate amber
        conda activate PROTACable
        grep "^HETATM" ${e3%.*}_prep1.pdb > ${e3%.*}_ligand0.pdb
        grep "^CONECT" ${e3%.*}_prep1.pdb >> ${e3%.*}_ligand0.pdb
        grep "^END" ${e3%.*}_prep1.pdb >> ${e3%.*}_ligand0.pdb

        sed -i '/HETATM/d' ${e3%.*}_prep1.pdb
        sed -i '/CONECT/d' ${e3%.*}_prep1.pdb
        sed -i '/END/d' ${e3%.*}_prep1.pdb

        pdbfixer ${e3%.*}_prep1.pdb --replace-nonstandard --add-residues --output="fixed_${e3%.*}_prep1.pdb"
        sed -i '/END/d' fixed_${e3%.*}_prep1.pdb
        sed -i '/TER.*UNL\|UNL.*TER/d' fixed_${e3%.*}_prep1.pdb
        cat ${e3%.*}_ligand0.pdb >> fixed_${e3%.*}_prep1.pdb

        conda deactivate PROTACable
        conda activate amber

        pdb4amber -i fixed_${e3%.*}_prep1.pdb -o ${e3%.*}_r_prep1.pdb
        cp $e3 ${e3%.*}_prep1.pdb
	

        #E3
        ###
        for i in {1..9}; do for atm in $o_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i}  UNL/}" ${e3%.*}_prep1.pdb;done;done
        for i in {10..80}; do for atm in $o_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i} UNL/}" ${e3%.*}_prep1.pdb;done;done

        for i in {1..9}; do for atm in $t_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i}  UNL/}" ${e3%.*}_prep1.pdb;done;done
        for i in {10..80}; do for atm in $t_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i} UNL/}" ${e3%.*}_prep1.pdb;done;done

        sed -i 's/ H  1 UNL/  H   UNL/g' ${e3%.*}_prep1.pdb; sed -i 's/ H  2 UNL/  H   UNL/g' ${e3%.*}_prep1.pdb; sed -i 's/ H  3 UNL/  H   UNL/g' ${e3%.*}_prep1.pdb;

        for hyd in $t_Hatoms; do sed -i "s/ $hyd  UNL/ H   UNL/g" ${e3%.*}_prep1.pdb;done
        for hyd in $th_Hatoms; do sed -i "s/ $hyd UNL/ H   UNL/g" ${e3%.*}_prep1.pdb;done
        for hyd in $f_Hatoms; do sed -i "s/ $hyd UNL/  H   UNL/g" ${e3%.*}_prep1.pdb;done

        for i in {1..9}; do sed -i "$((i-1)),/ H   UNL/{s/ H   UNL/ H${i}  UNL/}" ${e3%.*}_prep1.pdb;done
        for i in {10..80}; do sed -i "$((i-1)),/ H   UNL/{s/ H   UNL/ H${i} UNL/}" ${e3%.*}_prep1.pdb;done
	
        cp fixed_${e3%.*}_prep1.pdb ${e3%.*}_prep1_special.pdb
        ## spesh case
        rm ${e3%.*}_r_special.pdb ${e3%.*}_r_special2.pdb

        grep "^HETATM" ${e3%.*}_prep1_special.pdb > ${e3%.*}_prep1_special_lig.pdb
        grep "^CONECT" ${e3%.*}_prep1_special.pdb >> ${e3%.*}_prep1_special_lig.pdb
        grep "^END" ${e3%.*}_prep1_special.pdb >> ${e3%.*}_prep1_special_lig.pdb

        sed -i '/HETATM/d' ${e3%.*}_prep1_special.pdb
        sed -i '/CONECT/d' ${e3%.*}_prep1_special.pdb
        sed -i '/END/d' ${e3%.*}_prep1_special.pdb


        pdb4amber -i ${e3%.*}_prep1_special.pdb -o ${e3%.*}_r_special.pdb --nohyd
        sed -i '/HETATM/d' ${e3%.*}_r_special.pdb
        sed -i '/CONECT/d' ${e3%.*}_r_special.pdb
        sed -i '/END/d' ${e3%.*}_r_special.pdb

        cat ${e3%.*}_prep1_special_lig.pdb >> ${e3%.*}_r_special.pdb

        pdb4amber -i ${e3%.*}_r_special.pdb -o ${e3%.*}_r_special2.pdb
fi


echo "renaming HETATMs"

o_atoms="C O N F S P"
t_atoms="Cl Br"
o_Hatoms="H"
t_Hatoms="HN HO HS HP"
th_Hatoms="HXT HN1 HN2 HN3 HO1 1HO HO2 2HO HO3 3HO HS1 1HS HS2 2SH HS3 3SH HP1 1HP HP2 2HP HP3 3HP"
f_Hatoms="1HXT 2HXT 3HXT"

#POI
####
for i in {1..9}; do for atm in $o_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i}  UNL/}" ${poi%.*}_prep1.pdb;done;done
for i in {10..80}; do for atm in $o_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i} UNL/}" ${poi%.*}_prep1.pdb;done;done

for i in {1..9}; do for atm in $t_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i}  UNL/}" ${poi%.*}_prep1.pdb;done;done
for i in {10..80}; do for atm in $t_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i} UNL/}" ${poi%.*}_prep1.pdb;done;done

sed -i 's/ H  1 UNL/  H   UNL/g' ${poi%.*}_prep1.pdb; sed -i 's/ H  2 UNL/  H   UNL/g' ${poi%.*}_prep1.pdb; sed -i 's/ H  3 UNL/  H   UNL/g' ${poi%.*}_prep1.pdb;

for hyd in $t_Hatoms; do sed -i "s/ $hyd  UNL/ H   UNL/g" ${poi%.*}_prep1.pdb;done
for hyd in $th_Hatoms; do sed -i "s/ $hyd UNL/ H   UNL/g" ${poi%.*}_prep1.pdb;done
for hyd in $f_Hatoms; do sed -i "s/ $hyd UNL/  H   UNL/g" ${poi%.*}_prep1.pdb;done


for i in {1..9}; do sed -i "$((i-1)),/ H   UNL/{s/ H   UNL/ H${i}  UNL/}" ${poi%.*}_prep1.pdb;done
for i in {10..80}; do sed -i "$((i-1)),/ H   UNL/{s/ H   UNL/ H${i} UNL/}" ${poi%.*}_prep1.pdb;done

#E3
###
for i in {1..9}; do for atm in $o_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i}  UNL/}" ${e3%.*}_prep1.pdb;done;done
for i in {10..80}; do for atm in $o_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i} UNL/}" ${e3%.*}_prep1.pdb;done;done

for i in {1..9}; do for atm in $t_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i}  UNL/}" ${e3%.*}_prep1.pdb;done;done
for i in {10..80}; do for atm in $t_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i} UNL/}" ${e3%.*}_prep1.pdb;done;done

sed -i 's/ H  1 UNL/  H   UNL/g' ${e3%.*}_prep1.pdb; sed -i 's/ H  2 UNL/  H   UNL/g' ${e3%.*}_prep1.pdb; sed -i 's/ H  3 UNL/  H   UNL/g' ${e3%.*}_prep1.pdb;

for hyd in $t_Hatoms; do sed -i "s/ $hyd  UNL/ H   UNL/g" ${e3%.*}_prep1.pdb;done
for hyd in $th_Hatoms; do sed -i "s/ $hyd UNL/ H   UNL/g" ${e3%.*}_prep1.pdb;done
for hyd in $f_Hatoms; do sed -i "s/ $hyd UNL/  H   UNL/g" ${e3%.*}_prep1.pdb;done

for i in {1..9}; do sed -i "$((i-1)),/ H   UNL/{s/ H   UNL/ H${i}  UNL/}" ${e3%.*}_prep1.pdb;done
for i in {10..80}; do sed -i "$((i-1)),/ H   UNL/{s/ H   UNL/ H${i} UNL/}" ${e3%.*}_prep1.pdb;done


echo "initiating docking for complexes: $poi and $e3"

grep "^HETATM" ${poi%.*}_prep1.pdb > ${poi%.*}_ligand.pdb
#grep "^CONECT" $poi >> ${poi%.*}_ligand.pdb
grep "^HETATM" ${e3%.*}_prep1.pdb > ${e3%.*}_ligand.pdb
#grep "^CONECT" $e3 >> ${e3%.*}_ligand.pdb

sed -i '/UNL/!d' ${poi%.*}_ligand.pdb

sed -i '/UNL/!d' ${e3%.*}_ligand.pdb


sed -i '/ACE/d' ${poi%.*}_r_prep1.pdb
sed -i '/ACE/d' ${e3%.*}_r_prep1.pdb
sed -i '/NME/d' ${poi%.*}_r_prep1.pdb
sed -i '/NME/d' ${e3%.*}_r_prep1.pdb

sed -i '/HETATM/s/CL/Cl/g' ${poi%.*}_r_prep1.pdb
sed -i '/HETATM/s/CL/Cl/g' ${e3%.*}_r_prep1.pdb
sed -i '/HETATM/s/BR/Br/g' ${poi%.*}_r_prep1.pdb
sed -i '/HETATM/s/BR/Br/g' ${e3%.*}_r_prep1.pdb

############################################################################################################################################################
############################################################################################################################################################
sed -i '/^ATOM/s/CSD/CYS/g' ${poi%.*}_r_prep1.pdb                 #3-SULFINOALANINE
sed -i '/^ATOM/s/HYP/PRO/g' ${poi%.*}_r_prep1.pdb                 #4-HYDROXYPROLINE
sed -i '/^ATOM/s/BMT/THR/g' ${poi%.*}_r_prep1.pdb                 #4-METHYL-4-[(E)-2-BUTENYL]-4,N-METHYL-THREONINE
sed -i '/^ATOM/s/5HP/GLU/g' ${poi%.*}_r_prep1.pdb                 #5-HYDROXYPROLINE
sed -i '/^ATOM/s/ABA/ALA/g' ${poi%.*}_r_prep1.pdb                 #ALPHA-AMINOBUTYRIC_ACID
sed -i '/^ATOM/s/AIB/ALA/g' ${poi%.*}_r_prep1.pdb                 #ALPHA-AMINOISOBUTYRIC_ACID
sed -i '/^ATOM/s/CSW/CYS/g' ${poi%.*}_r_prep1.pdb                 #CYSTEINE-S-DIOXIDE
sed -i '/^ATOM/s/OCS/CYS/g' ${poi%.*}_r_prep1.pdb                 #CYSTEINESULFONIC_ACID
sed -i '/^ATOM/s/DAL/ALA/g' ${poi%.*}_r_prep1.pdb                 #D-ALANINE
sed -i '/^ATOM/s/DAR/ARG/g' ${poi%.*}_r_prep1.pdb                 #D-ARGININE
sed -i '/^ATOM/s/DSG/ASN/g' ${poi%.*}_r_prep1.pdb                 #D-ASPARAGINE
sed -i '/^ATOM/s/DSP/ASP/g' ${poi%.*}_r_prep1.pdb                 #D-ASPARTATE
sed -i '/^ATOM/s/DCY/CYS/g' ${poi%.*}_r_prep1.pdb                 #D-CYSTEINE
sed -i '/^ATOM/s/CRO/CRO/g' ${poi%.*}_r_prep1.pdb                 #DECARBOXY(PARAHYDROXYBENZYLIDENE-IMIDAZOLIDINONE)THREONINE
sed -i '/^ATOM/s/DGL/GLU/g' ${poi%.*}_r_prep1.pdb                 #D-GLUTAMATE
sed -i '/^ATOM/s/DGN/GLN/g' ${poi%.*}_r_prep1.pdb                 #D-GLUTAMINE
sed -i '/^ATOM/s/DHI/HIS/g' ${poi%.*}_r_prep1.pdb                 #D-HISTIDINE
sed -i '/^ATOM/s/DIL/ILE/g' ${poi%.*}_r_prep1.pdb                 #D-ISOLEUCINE
sed -i '/^ATOM/s/DIV/VAL/g' ${poi%.*}_r_prep1.pdb                 #D-ISOVALINE
sed -i '/^ATOM/s/DLE/LEU/g' ${poi%.*}_r_prep1.pdb                 #D-LEUCINE
sed -i '/^ATOM/s/DLY/LYS/g' ${poi%.*}_r_prep1.pdb                 #D-LYSINE
sed -i '/^ATOM/s/DPN/PHE/g' ${poi%.*}_r_prep1.pdb                 #D-PHENYLALANINE
sed -i '/^ATOM/s/DPR/PRO/g' ${poi%.*}_r_prep1.pdb                 #D-PROLINE
sed -i '/^ATOM/s/DSN/SER/g' ${poi%.*}_r_prep1.pdb                 #D-SERINE
sed -i '/^ATOM/s/DTH/THR/g' ${poi%.*}_r_prep1.pdb                 #D-THREONINE
sed -i '/^ATOM/s/DTR/DTR/g' ${poi%.*}_r_prep1.pdb                 #D-TRYPTOPHANE
sed -i '/^ATOM/s/DTY/TYR/g' ${poi%.*}_r_prep1.pdb                 #D-TYROSINE
sed -i '/^ATOM/s/DVA/VAL/g' ${poi%.*}_r_prep1.pdb                 #D-VALINE
sed -i '/^ATOM/s/CGU/GLU/g' ${poi%.*}_r_prep1.pdb                 #GAMMA-CARBOXY-GLUTAMIC_ACID
sed -i '/^ATOM/s/KCX/LYS/g' ${poi%.*}_r_prep1.pdb                 #LYSINE_NZ-CARBOXYLIC_ACID
sed -i '/^ATOM/s/LLP/LYS/g' ${poi%.*}_r_prep1.pdb                 #LYSINE-PYRIDOXAL-5'-PHOSPHATE
sed -i '/^ATOM/s/CXM/MET/g' ${poi%.*}_r_prep1.pdb                 #N-CARBOXYMETHIONINE
sed -i '/^ATOM/s/FME/MET/g' ${poi%.*}_r_prep1.pdb                 #N-FORMYLMETHIONINE
sed -i '/^ATOM/s/MLE/LEU/g' ${poi%.*}_r_prep1.pdb                 #N-METHYLLEUCINE
sed -i '/^ATOM/s/MVA/VAL/g' ${poi%.*}_r_prep1.pdb                 #N-METHYLVALINE
sed -i '/^ATOM/s/NLE/LEU/g' ${poi%.*}_r_prep1.pdb                 #NORLEUCINE
sed -i '/^ATOM/s/PTR/TYR/g' ${poi%.*}_r_prep1.pdb                 #O-PHOSPHOTYROSINE
sed -i '/^ATOM/s/ORN/ALA/g' ${poi%.*}_r_prep1.pdb                 #ORNITHINE
sed -i '/^ATOM/s/SEP/SER/g' ${poi%.*}_r_prep1.pdb                 #PHOSPHOSERINE
sed -i '/^ATOM/s/TPO/THR/g' ${poi%.*}_r_prep1.pdb                 #PHOSPHOTHREONINE
sed -i '/^ATOM/s/PCA/GLU/g' ${poi%.*}_r_prep1.pdb                 #PYROGLUTAMIC_ACID
sed -i '/^ATOM/s/SAR/GLY/g' ${poi%.*}_r_prep1.pdb                 #SARCOSINE
sed -i '/^ATOM/s/CEA/CYS/g' ${poi%.*}_r_prep1.pdb                 #S-HYDROXY-CYSTEINE
sed -i '/^ATOM/s/CSO/CYS/g' ${poi%.*}_r_prep1.pdb                 #S-HYDROXYCYSTEINE
sed -i '/^ATOM/s/CSS/CYS/g' ${poi%.*}_r_prep1.pdb                 #S-MERCAPTOCYSTEINE
sed -i '/^ATOM/s/CSX/CYS/g' ${poi%.*}_r_prep1.pdb                 #S-OXY_CYSTEINE
sed -i '/^ATOM/s/CME/CYS/g' ${poi%.*}_r_prep1.pdb                 #S,S-(2-HYDROXYETHYL)THIOCYSTEINE
sed -i '/^ATOM/s/TYS/TYR/g' ${poi%.*}_r_prep1.pdb                 #SULFONATED_TYROSINE
sed -i '/^ATOM/s/TPQ/PHE/g' ${poi%.*}_r_prep1.pdb                 #TOPO-QUINONE
sed -i '/^ATOM/s/STY/TYR/g' ${poi%.*}_r_prep1.pdb                 #TYROSINE-O-SULPHONIC_ACID
sed -i '/^ATOM/s/CCS/CYS/g' ${poi%.*}_r_prep1.pdb                 #https://github.com/EvanBaugh/miscellaneous_scripts/blob/master/process_pdb.py
sed -i '/^ATOM/s/CALA/ALA/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CARG/ARG/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CASN/ASN/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CASP/ASP/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CCYS/CYS/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CCYX/CYX/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CGLN/GLN/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CGLU/GLU/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CGLY/GLY/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CHID/HID/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CHIE/HIE/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CHIP/HIP/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CHYP/HYP/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CILE/ILE/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CLEU/LEU/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CLYS/LYS/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CMET/MET/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CPHE/PHE/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CPRO/PRO/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CSER/SER/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CTHR/THR/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CTRP/TRP/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CTYR/TYR/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CVAL/VAL/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NALA/ALA/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NARG/ARG/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NASN/ASN/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NASP/ASP/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NCYS/CYS/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NCYX/CYX/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NGLN/GLN/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NGLU/GLU/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NGLY/GLY/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NHID/HID/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NHIE/HIE/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NHIP/HIP/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NILE/ILE/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NLEU/LEU/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NLYS/LYS/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NMET/MET/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NPHE/PHE/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NPRO/PRO/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NSER/SER/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NTHR/THR/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NTRP/TRP/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NTYR/TYR/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NVAL/VAL/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/DAS/ASP/g' ${poi%.*}_r_prep1.pdb                 #GUESS https://github.com/EvanBaugh/miscellaneous_scripts/blob/master/process_pdb.py
sed -i '/^ATOM/s/CAF/CYS/g' ${poi%.*}_r_prep1.pdb                 #GUESS https://github.com/EvanBaugh/miscellaneous_scripts/blob/master/process_pdb.py ## S-DIMETHYLARSINOYL-CYSTEINE
sed -i '/^ATOM/s/CAS/CYS/g' ${poi%.*}_r_prep1.pdb                 #GUESS https://github.com/EvanBaugh/miscellaneous_scripts/blob/master/process_pdb.py ## S-(DIMETHYLARSENIC)CYSTEINE
sed -i '/^ATOM/s/AIB/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/ALA/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  ALA
sed -i '/^ATOM/s/ALM/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/AYA/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/BNN/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/CHG/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/CSD/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/DAL/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/DHA/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/DNP/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/FLA/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/HAC/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/PRR/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/MAA/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/TIH/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/TPQ/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/0CS/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    0CS ALA  3-[(S)-HYDROPEROXYSULFINYL]-L-ALANINE
sed -i '/^ATOM/s/2BU/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    2BU ADE
sed -i '/^ATOM/s/2OP/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    2OP (2S  2-HYDROXYPROPANAL
sed -i '/^ATOM/s/4F3/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    4F3 ALA  CYCLIZED
sed -i '/^ATOM/s/AA4/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    AA4 ALA  2-AMINO-5-HYDROXYPENTANOIC ACID
sed -i '/^ATOM/s/ABA/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    ABA ALA  ALPHA-AMINOBUTYRIC ACID
sed -i '/^ATOM/s/AHO/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    AHO ALA  N-ACETYL-N-HYDROXY-L-ORNITHINE
sed -i '/^ATOM/s/AHP/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    AHP ALA  2-AMINO-HEPTANOIC ACID
sed -i '/^ATOM/s/AIB/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    AIB ALA  ALPHA-AMINOISOBUTYRIC ACID
sed -i '/^ATOM/s/ALA/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    ALA ALA
sed -i '/^ATOM/s/ALC/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    ALC ALA  2-AMINO-3-CYCLOHEXYL-PROPIONIC ACID
sed -i '/^ATOM/s/ALM/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    ALM ALA  1-METHYL-ALANINAL
sed -i '/^ATOM/s/ALN/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    ALN ALA  NAPHTHALEN-2-YL-3-ALANINE
sed -i '/^ATOM/s/ALS/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    ALS ALA  2-AMINO-3-OXO-4-SULFO-BUTYRIC ACID
sed -i '/^ATOM/s/ALT/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    ALT ALA  THIOALANINE
sed -i '/^ATOM/s/AP7/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    AP7 ADE
sed -i '/^ATOM/s/APH/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    APH ALA  P-AMIDINOPHENYL-3-ALANINE
sed -i '/^ATOM/s/AYA/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    AYA ALA  N-ACETYLALANINE
sed -i '/^ATOM/s/AYG/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    AYG ALA
sed -i '/^ATOM/s/B2A/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    B2A ALA  ALANINE BORONIC ACID
sed -i '/^ATOM/s/B3A/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    B3A ALA  (3S)-3-AMINOBUTANOIC ACID
sed -i '/^ATOM/s/BAL/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    BAL ALA  BETA-ALANINE
sed -i '/^ATOM/s/BNN/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    BNN ALA  ACETYL-P-AMIDINOPHENYLALANINE
sed -i '/^ATOM/s/C12/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    C12 ALA
sed -i '/^ATOM/s/C99/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    C99 ALA
sed -i '/^ATOM/s/CAB/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CAB ALA  4-CARBOXY-4-AMINOBUTANAL
sed -i '/^ATOM/s/CH6/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CH6 ALA
sed -i '/^ATOM/s/CH7/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CH7 ALA
sed -i '/^ATOM/s/CLB/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CLB ALA
sed -i '/^ATOM/s/CLD/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CLD ALA
sed -i '/^ATOM/s/CLV/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CLV ALA
sed -i '/^ATOM/s/CQR/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CQR ALA
sed -i '/^ATOM/s/CR2/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CR2 ALA  POST-TRANSLATIONAL MODIFICATION
sed -i '/^ATOM/s/CR5/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CR5 ALA
sed -i '/^ATOM/s/CR7/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CR7 ALA
sed -i '/^ATOM/s/CR8/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CR8 ALA
sed -i '/^ATOM/s/CRK/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CRK ALA
sed -i '/^ATOM/s/CRW/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CRW ALA
sed -i '/^ATOM/s/CRX/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CRX ALA
sed -i '/^ATOM/s/CSI/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CSI ALA
sed -i '/^ATOM/s/CSY/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CSY ALA  MODIFIED TYROSINE COMPLEX
sed -i '/^ATOM/s/CWR/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CWR ALA
sed -i '/^ATOM/s/DAB/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    DAB ALA  24-DIAMINOBUTYRIC ACID
sed -i '/^ATOM/s/DAL/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    DAL ALA  D-ALANINE
sed -i '/^ATOM/s/DAM/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    DAM ALA  N-METHYL-ALPHA-BETA-DEHYDROALANINE
sed -i '/^ATOM/s/DBU/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    DBU ALA  (2E)-2-AMINOBUT-2-ENOIC ACID
sed -i '/^ATOM/s/DBZ/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    DBZ ALA  3-(BENZOYLAMINO)-L-ALANINE
sed -i '/^ATOM/s/DHA/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    DHA ALA  2-AMINO-ACRYLIC ACID
sed -i '/^ATOM/s/DPP/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    DPP ALA  DIAMMINOPROPANOIC ACID
sed -i '/^ATOM/s/FGL/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    FGL ALA  2-AMINOPROPANEDIOIC ACID
sed -i '/^ATOM/s/HHK/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    HHK ALA  (2S)-28-DIAMINOOCTANOIC ACID
sed -i '/^ATOM/s/HMF/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    HMF ALA  2-AMINO-4-PHENYL-BUTYRIC ACID
sed -i '/^ATOM/s/IAM/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    IAM ALA  4-[(ISOPROPYLAMINO)METHYL]PHENYLALANINE
sed -i '/^ATOM/s/IGL/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    IGL ALA  ALPHA-AMINO-2-INDANACETIC ACID
sed -i '/^ATOM/s/KYN/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    KYN ALA  KYNURENINE
sed -i '/^ATOM/s/LAL/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    LAL ALA  NN-DIMETHYL-L-ALANINE
sed -i '/^ATOM/s/MAA/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    MAA ALA  N-METHYLALANINE
sed -i '/^ATOM/s/MDO/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    MDO ALA
sed -i '/^ATOM/s/MFC/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    MFC ALA  CYCLIZED
sed -i '/^ATOM/s/NAL/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    NAL ALA  BETA-(2-NAPHTHYL)-ALANINE
sed -i '/^ATOM/s/NAM/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    NAM ALA  NAM NAPTHYLAMINOALANINE
sed -i '/^ATOM/s/NCB/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    NCB ALA  CHEMICAL MODIFICATION
sed -i '/^ATOM/s/NRQ/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    NRQ ALA
sed -i '/^ATOM/s/NYC/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    NYC ALA
sed -i '/^ATOM/s/ORN/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    ORN ALA  ORNITHINE
sed -i '/^ATOM/s/PIA/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    PIA ALA  FUSION OF ALA 65 TYR 66 GLY 67
sed -i '/^ATOM/s/PRR/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    PRR ALA  3-(METHYL-PYRIDINIUM)ALANINE
sed -i '/^ATOM/s/PYA/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    PYA ALA  3-(110-PHENANTHROL-2-YL)-L-ALANINE
sed -i '/^ATOM/s/PYC/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    PYC ALA  PYRROLE-2-CARBOXYLATE
sed -i '/^ATOM/s/PYT/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    PYT ALA  MODIFIED ALANINE
sed -i '/^ATOM/s/RC7/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    RC7 ALA
sed -i '/^ATOM/s/SEC/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    SEC ALA  2-AMINO-3-SELENINO-PROPIONIC ACID
sed -i '/^ATOM/s/SIC/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    SIC ALA
sed -i '/^ATOM/s/SUI/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    SUI ALA
sed -i '/^ATOM/s/TIH/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    TIH ALA  BETA(2-THIENYL)ALANINE
sed -i '/^ATOM/s/TPQ/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    TPQ ALA  245-TRIHYDROXYPHENYLALANINE
sed -i '/^ATOM/s/UMA/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    UMA ALA
sed -i '/^ATOM/s/X9Q/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    X9Q ALA
sed -i '/^ATOM/s/XXY/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    XXY ALA
sed -i '/^ATOM/s/XYG/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    XYG ALA
sed -i '/^ATOM/s/BCS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/BUC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/C5C/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/C6C/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/CCS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/CEA/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/CME/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/CSO/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/CSP/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/CSS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/CSX/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/CSW/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/CY1/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/CY3/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/CYG/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/CYM/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/CYS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  CYS
sed -i '/^ATOM/s/CYQ/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/DCY/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/EFC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/OCS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/PEC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/PR3/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/SCH/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/SCS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/SCY/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/SHC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/SMC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/SOC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/5CS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    5CS CYS
sed -i '/^ATOM/s/AGT/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    AGT CYS  AGMATINE-CYSTEINE ADDUCT
sed -i '/^ATOM/s/BBC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    BBC CYS
sed -i '/^ATOM/s/BCS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    BCS CYS  BENZYLCYSTEINE
sed -i '/^ATOM/s/BCX/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    BCX CYS  BETA-3-CYSTEINE
sed -i '/^ATOM/s/BPE/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    BPE CYS
sed -i '/^ATOM/s/BUC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    BUC CYS  SS-BUTYLTHIOCYSTEINE
sed -i '/^ATOM/s/C3Y/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    C3Y CYS  MODIFIED CYSTEINE
sed -i '/^ATOM/s/C5C/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    C5C CYS  S-CYCLOPENTYL THIOCYSTEINE
sed -i '/^ATOM/s/C6C/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    C6C CYS  S-CYCLOHEXYL THIOCYSTEINE
sed -i '/^ATOM/s/CAF/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CAF CYS  S-DIMETHYLARSINOYL-CYSTEINE
sed -i '/^ATOM/s/CAS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CAS CYS  S-(DIMETHYLARSENIC)CYSTEINE
sed -i '/^ATOM/s/CCS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CCS CYS  CARBOXYMETHYLATED CYSTEINE
sed -i '/^ATOM/s/CME/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CME CYS  MODIFIED CYSTEINE
sed -i '/^ATOM/s/CML/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CML CYS
sed -i '/^ATOM/s/CMT/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CMT CYS  O-METHYLCYSTEINE
sed -i '/^ATOM/s/CS1/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CS1 CYS  S-(2-ANILINYL-SULFANYL)-CYSTEINE
sed -i '/^ATOM/s/CS3/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CS3 CYS
sed -i '/^ATOM/s/CS4/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CS4 CYS
sed -i '/^ATOM/s/CSA/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CSA CYS  S-ACETONYLCYSTEIN
sed -i '/^ATOM/s/CSB/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CSB CYS  CYS BOUND TO LEAD ION
sed -i '/^ATOM/s/CSD/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CSD CYS  3-SULFINOALANINE
sed -i '/^ATOM/s/CSE/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CSE CYS  SELENOCYSTEINE
sed -i '/^ATOM/s/CSO/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CSO CYS  INE S-HYDROXYCYSTEINE
sed -i '/^ATOM/s/CSR/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CSR CYS  S-ARSONOCYSTEINE
sed -i '/^ATOM/s/CSS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CSS CYS  13-THIAZOLE-4-CARBOXYLIC ACID
sed -i '/^ATOM/s/CSU/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CSU CYS  CYSTEINE-S-SULFONIC ACID
sed -i '/^ATOM/s/CSW/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CSW CYS  CYSTEINE-S-DIOXIDE
sed -i '/^ATOM/s/CSX/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CSX CYS  OXOCYSTEINE
sed -i '/^ATOM/s/CSZ/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CSZ CYS  S-SELANYL CYSTEINE
sed -i '/^ATOM/s/CY0/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CY0 CYS  MODIFIED CYSTEINE
sed -i '/^ATOM/s/CY1/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CY1 CYS  ACETAMIDOMETHYLCYSTEINE
sed -i '/^ATOM/s/CY3/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CY3 CYS  2-AMINO-3-MERCAPTO-PROPIONAMIDE
sed -i '/^ATOM/s/CY4/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CY4 CYS  S-BUTYRYL-CYSTEIN
sed -i '/^ATOM/s/CY7/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CY7 CYS  MODIFIED CYSTEINE
sed -i '/^ATOM/s/CYD/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CYD CYS
sed -i '/^ATOM/s/CYF/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CYF CYS  FLUORESCEIN LABELLED CYS380 (P14)
sed -i '/^ATOM/s/CYG/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CYG CYS
sed -i '/^ATOM/s/CYQ/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CYQ CYS
sed -i '/^ATOM/s/CYR/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CYR CYS
sed -i '/^ATOM/s/CYS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CYS CYS
sed -i '/^ATOM/s/CZ2/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CZ2 CYS  S-(DIHYDROXYARSINO)CYSTEINE
sed -i '/^ATOM/s/CZZ/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CZZ CYS  THIARSAHYDROXY-CYSTEINE
sed -i '/^ATOM/s/DCY/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    DCY CYS  D-CYSTEINE
sed -i '/^ATOM/s/DYS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    DYS CYS
sed -i '/^ATOM/s/EFC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    EFC CYS  SS-(2-FLUOROETHYL)THIOCYSTEINE
sed -i '/^ATOM/s/FOE/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    FOE CYS
sed -i '/^ATOM/s/GT9/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    GT9 CYS  SG ALKYLATED
sed -i '/^ATOM/s/GYC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    GYC CYS
sed -i '/^ATOM/s/HTI/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    HTI CYS
sed -i '/^ATOM/s/KOR/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    KOR CYS  MODIFIED CYSTEINE
sed -i '/^ATOM/s/M0H/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    M0H CYS  S-(HYDROXYMETHYL)-L-CYSTEINE
sed -i '/^ATOM/s/MCS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    MCS CYS  MALONYLCYSTEINE
sed -i '/^ATOM/s/NPH/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    NPH CYS
sed -i '/^ATOM/s/NYS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    NYS CYS
sed -i '/^ATOM/s/OCS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    OCS CYS  CYSTEINE SULFONIC ACID
sed -i '/^ATOM/s/OCY/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    OCY CYS  HYDROXYETHYLCYSTEINE
sed -i '/^ATOM/s/P1L/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    P1L CYS  S-PALMITOYL CYSTEINE
sed -i '/^ATOM/s/PBB/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    PBB CYS  S-(4-BROMOBENZYL)CYSTEINE
sed -i '/^ATOM/s/PEC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    PEC CYS  SS-PENTYLTHIOCYSTEINE
sed -i '/^ATOM/s/PR3/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    PR3 CYS  INE DTT-CYSTEINE
sed -i '/^ATOM/s/PYX/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    PYX CYS  S-[S-THIOPYRIDOXAMINYL]CYSTEINE
sed -i '/^ATOM/s/R1A/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    R1A CYS
sed -i '/^ATOM/s/R1B/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    R1B CYS
sed -i '/^ATOM/s/R1F/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    R1F CYS
sed -i '/^ATOM/s/R7A/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    R7A CYS
sed -i '/^ATOM/s/RCY/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    RCY CYS
sed -i '/^ATOM/s/SAH/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    SAH CYS  S-ADENOSYL-L-HOMOCYSTEINE
sed -i '/^ATOM/s/SC2/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    SC2 CYS  N-ACETYL-L-CYSTEINE
sed -i '/^ATOM/s/SCH/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    SCH CYS  S-METHYL THIOCYSTEINE GROUP
sed -i '/^ATOM/s/SCS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    SCS CYS  MODIFIED CYSTEINE
sed -i '/^ATOM/s/SCY/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    SCY CYS  CETYLATED CYSTEINE
sed -i '/^ATOM/s/SHC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    SHC CYS  S-HEXYLCYSTEINE
sed -i '/^ATOM/s/SMC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    SMC CYS  POST-TRANSLATIONAL MODIFICATION
sed -i '/^ATOM/s/SNC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    SNC CYS  S-NITROSO CYSTEINE
sed -i '/^ATOM/s/SOC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    SOC CYS  DIOXYSELENOCYSTEINE
sed -i '/^ATOM/s/TEE/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    TEE CYS  POST-TRANSLATIONAL MODIFICATION
sed -i '/^ATOM/s/TNB/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    TNB CYS  S-(236-TRINITROPHENYL)CYSTEINE
sed -i '/^ATOM/s/TYX/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    TYX CYS  S-(2-ANILINO-2-OXOETHYL)-L-CYSTEINE
sed -i '/^ATOM/s/YCM/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    YCM CYS  S-(2-AMINO-2-OXOETHYL)-L-CYSTEINE
sed -i '/^ATOM/s/2AS/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
sed -i '/^ATOM/s/ASA/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
sed -i '/^ATOM/s/ASB/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
sed -i '/^ATOM/s/ASK/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
sed -i '/^ATOM/s/ASL/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
sed -i '/^ATOM/s/ASP/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                  ASP
sed -i '/^ATOM/s/ASQ/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
sed -i '/^ATOM/s/BHD/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
sed -i '/^ATOM/s/DAS/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
sed -i '/^ATOM/s/DSP/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
sed -i '/^ATOM/s/3MD/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    3MD ASP  2S3S-3-METHYLASPARTIC ACID
sed -i '/^ATOM/s/A0A/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    A0A ASP  ASPARTYL-FORMYL MIXED ANHYDRIDE
sed -i '/^ATOM/s/ACB/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    ACB ASP  3-METHYL-ASPARTIC ACID
sed -i '/^ATOM/s/AKL/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    AKL ASP  3-AMINO-5-CHLORO-4-OXOPENTANOIC ACID
sed -i '/^ATOM/s/ASA/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    ASA ASP  ASPARTIC ALDEHYDE
sed -i '/^ATOM/s/ASB/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    ASB ASP  ASPARTIC ACID-4-CARBOXYETHYL ESTER
sed -i '/^ATOM/s/ASI/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    ASI ASP  L-ISO-ASPARTATE
sed -i '/^ATOM/s/ASK/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    ASK ASP  DEHYDROXYMETHYLASPARTIC ACID
sed -i '/^ATOM/s/ASL/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    ASL ASP  ASPARTIC ACID-4-CARBOXYETHYL ESTER
sed -i '/^ATOM/s/ASP/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    ASP ASP
sed -i '/^ATOM/s/B3D/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    B3D ASP  3-AMINOPENTANEDIOIC ACID
sed -i '/^ATOM/s/BFD/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    BFD ASP  ASPARTATE BERYLLIUM FLUORIDE
sed -i '/^ATOM/s/BHD/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    BHD ASP  BETA-HYDROXYASPARTIC ACID
sed -i '/^ATOM/s/DAS/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    DAS ASP  D-ASPARTIC ACID
sed -i '/^ATOM/s/DMK/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    DMK ASP  DIMETHYL ASPARTIC ACID
sed -i '/^ATOM/s/IAS/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    IAS ASP  ASPARTYL GROUP
sed -i '/^ATOM/s/OHS/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    OHS ASP  O-(CARBOXYSULFANYL)-4-OXO-L-HOMOSERINE
sed -i '/^ATOM/s/OXX/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    OXX ASP  OXALYL-ASPARTYL ANHYDRIDE
sed -i '/^ATOM/s/PHD/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    PHD ASP  2-AMINO-4-OXO-4-PHOSPHONOOXY-BUTYRIC ACID
sed -i '/^ATOM/s/SNN/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    SNN ASP  POST-TRANSLATIONAL MODIFICATION
sed -i '/^ATOM/s/5HP/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
sed -i '/^ATOM/s/CGU/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
sed -i '/^ATOM/s/DGL/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
sed -i '/^ATOM/s/GGL/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
sed -i '/^ATOM/s/GLU/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                  GLU
sed -i '/^ATOM/s/GMA/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
sed -i '/^ATOM/s/PCA/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
sed -i '/^ATOM/s/AB7/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                    AB7 GLU  ALPHA-AMINOBUTYRIC ACID
sed -i '/^ATOM/s/AR4/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                    AR4 GLU
sed -i '/^ATOM/s/B3E/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                    B3E GLU  (3S)-3-AMINOHEXANEDIOIC ACID
sed -i '/^ATOM/s/CGU/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                    CGU GLU  CARBOXYLATION OF THE CG ATOM
sed -i '/^ATOM/s/DGL/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                    DGL GLU  D-GLU
sed -i '/^ATOM/s/GLU/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                    GLU GLU
sed -i '/^ATOM/s/GMA/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                    GMA GLU  1-AMIDO-GLUTAMIC ACID
sed -i '/^ATOM/s/ILG/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                    ILG GLU  GLU LINKED TO NEXT RESIDUE VIA CG
sed -i '/^ATOM/s/LME/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                    LME GLU  (3R)-3-METHYL-L-GLUTAMIC ACID
sed -i '/^ATOM/s/MEG/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                    MEG GLU  (2S3R)-3-METHYL-GLUTAMIC ACID
sed -i '/^ATOM/s/DAH/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
sed -i '/^ATOM/s/DPN/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
sed -i '/^ATOM/s/HPQ/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
sed -i '/^ATOM/s/PHE/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                  PHE
sed -i '/^ATOM/s/PHI/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
sed -i '/^ATOM/s/PHL/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
sed -i '/^ATOM/s/1PA/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    1PA PHE  PHENYLMETHYLACETIC ACID ALANINE
sed -i '/^ATOM/s/23F/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    23F PHE  (2Z)-2-AMINO-3-PHENYLACRYLIC ACID
sed -i '/^ATOM/s/4PH/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    4PH PHE  4-METHYL-L-PHENYLALANINE
sed -i '/^ATOM/s/B2F/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    B2F PHE  PHENYLALANINE BORONIC ACID
sed -i '/^ATOM/s/BIF/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    BIF PHE
sed -i '/^ATOM/s/CHS/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    CHS PHE  4-AMINO-5-CYCLOHEXYL-3-HYDROXY-PENTANOIC AC
sed -i '/^ATOM/s/DAH/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    DAH PHE  34-DIHYDROXYDAHNYLALANINE
sed -i '/^ATOM/s/DPH/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    DPH PHE  DEAMINO-METHYL-PHENYLALANINE
sed -i '/^ATOM/s/DPN/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    DPN PHE  D-CONFIGURATION
sed -i '/^ATOM/s/FCL/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    FCL PHE  3-CHLORO-L-PHENYLALANINE
sed -i '/^ATOM/s/FOG/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    FOG PHE  PHENYLALANINOYL-[1-HYDROXY]-2-PROPYLENE
sed -i '/^ATOM/s/FRF/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    FRF PHE  PHE FOLLOWED BY REDUCED PHE
sed -i '/^ATOM/s/HPE/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    HPE PHE  HOMOPHENYLALANINE
sed -i '/^ATOM/s/HPH/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    HPH PHE  PHENYLALANINOL GROUP
sed -i '/^ATOM/s/HPQ/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    HPQ PHE  HOMOPHENYLALANINYLMETHANE
sed -i '/^ATOM/s/MEA/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    MEA PHE  N-METHYLPHENYLALANINE
sed -i '/^ATOM/s/MTY/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    MTY PHE  3-HYDROXYPHENYLALANINE
sed -i '/^ATOM/s/NFA/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    NFA PHE  MODIFIED PHENYLALANINE
sed -i '/^ATOM/s/PBF/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    PBF PHE  PARA-(BENZOYL)-PHENYLALANINE
sed -i '/^ATOM/s/PCS/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    PCS PHE  PHENYLALANYLMETHYLCHLORIDE
sed -i '/^ATOM/s/PF5/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    PF5 PHE  23456-PENTAFLUORO-L-PHENYLALANINE
sed -i '/^ATOM/s/PFF/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    PFF PHE  4-FLUORO-L-PHENYLALANINE
sed -i '/^ATOM/s/PHA/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    PHA PHE  PHENYLALANINAL
sed -i '/^ATOM/s/PHE/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    PHE PHE
sed -i '/^ATOM/s/PHI/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    PHI PHE  IODO-PHENYLALANINE
sed -i '/^ATOM/s/PHL/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    PHL PHE  L-PHENYLALANINOL
sed -i '/^ATOM/s/PHM/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    PHM PHE  PHENYLALANYLMETHANE
sed -i '/^ATOM/s/PM3/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    PM3 PHE
sed -i '/^ATOM/s/PPN/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    PPN PHE  THE LIGAND IS A PARA-NITRO-PHENYLALANINE
sed -i '/^ATOM/s/PRQ/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    PRQ PHE  PHENYLALANINE
sed -i '/^ATOM/s/PSA/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    PSA PHE
sed -i '/^ATOM/s/SMF/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    SMF PHE  4-SULFOMETHYL-L-PHENYLALANINE
sed -i '/^ATOM/s/GL3/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
sed -i '/^ATOM/s/GLY/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                  GLY
sed -i '/^ATOM/s/GLZ/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
sed -i '/^ATOM/s/GSC/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
sed -i '/^ATOM/s/MPQ/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
sed -i '/^ATOM/s/MSA/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
sed -i '/^ATOM/s/NMC/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
sed -i '/^ATOM/s/SAR/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
sed -i '/^ATOM/s/ACY/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    ACY GLY  POST-TRANSLATIONAL MODIFICATION
sed -i '/^ATOM/s/CHG/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    CHG GLY  CYCLOHEXYL GLYCINE
sed -i '/^ATOM/s/CHP/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    CHP GLY  3-CHLORO-4-HYDROXYPHENYLGLYCINE
sed -i '/^ATOM/s/GHP/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    GHP GLY  4-HYDROXYPHENYLGLYCINE
sed -i '/^ATOM/s/GL3/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    GL3 GLY  POST-TRANSLATIONAL MODIFICATION
sed -i '/^ATOM/s/GLY/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    GLY GLY
sed -i '/^ATOM/s/GLZ/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    GLZ GLY  AMINO-ACETALDEHYDE
sed -i '/^ATOM/s/GYS/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    GYS GLY
sed -i '/^ATOM/s/IPG/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    IPG GLY  N-ISOPROPYL GLYCINE
sed -i '/^ATOM/s/MEU/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    MEU GLY  O-METHYL-GLYCINE
sed -i '/^ATOM/s/MPQ/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    MPQ GLY  N-METHYL-ALPHA-PHENYL-GLYCINE
sed -i '/^ATOM/s/MSA/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    MSA GLY  (2-S-METHYL) SARCOSINE
sed -i '/^ATOM/s/NMC/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    NMC GLY  N-CYCLOPROPYLMETHYL GLYCINE
sed -i '/^ATOM/s/PG9/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    PG9 GLY  D-PHENYLGLYCINE
sed -i '/^ATOM/s/SAR/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    SAR GLY  SARCOSINE
sed -i '/^ATOM/s/SHP/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    SHP GLY  (4-HYDROXYMALTOSEPHENYL)GLYCINE
sed -i '/^ATOM/s/TBG/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    TBG GLY  T-BUTYL GLYCINE
sed -i '/^ATOM/s/3AH/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
sed -i '/^ATOM/s/DHI/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
sed -i '/^ATOM/s/HIC/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
sed -i '/^ATOM/s/HIS/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                  HIS
sed -i '/^ATOM/s/MHS/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
sed -i '/^ATOM/s/NEM/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
sed -i '/^ATOM/s/NEP/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
sed -i '/^ATOM/s/HID/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                  single delta N protonation
sed -i '/^ATOM/s/HIE/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                  single epsilon N protonation
sed -i '/^ATOM/s/3AH/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                    3AH HIS
sed -i '/^ATOM/s/DDE/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                    DDE HIS
sed -i '/^ATOM/s/DHI/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                    DHI HIS  D-HISTIDINE
sed -i '/^ATOM/s/HIA/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                    HIA HIS  L-HISTIDINE AMIDE
sed -i '/^ATOM/s/HIC/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                    HIC HIS  4-METHYL-HISTIDINE
sed -i '/^ATOM/s/HIP/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                    HIP HIS  ND1-PHOSPHONOHISTIDINE...or commonly used doubly protonated state
sed -i '/^ATOM/s/HIQ/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                    HIQ HIS  MODIFIED HISTIDINE
sed -i '/^ATOM/s/HIS/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                    HIS HIS
sed -i '/^ATOM/s/HSO/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                    HSO HIS  HISTIDINOL
sed -i '/^ATOM/s/MHS/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                    MHS HIS  1-N-METHYLHISTIDINE
sed -i '/^ATOM/s/NEP/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                    NEP HIS  N1-PHOSPHONOHISTIDINE
sed -i '/^ATOM/s/NZH/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                    NZH HIS
sed -i '/^ATOM/s/OHI/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                    OHI HIS  3-(2-OXO-2H-IMIDAZOL-4-YL)-L-ALANINE
sed -i '/^ATOM/s/PSH/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                    PSH HIS  1-THIOPHOSPHONO-L-HISTIDINE
sed -i '/^ATOM/s/DIL/ILE/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ILE
sed -i '/^ATOM/s/IIL/ILE/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ILE
sed -i '/^ATOM/s/ILE/ILE/g' ${poi%.*}_r_prep1.pdb                 ###                  ILE
sed -i '/^ATOM/s/B2I/ILE/g' ${poi%.*}_r_prep1.pdb                 ###                    B2I ILE  ISOLEUCINE BORONIC ACID
sed -i '/^ATOM/s/DIL/ILE/g' ${poi%.*}_r_prep1.pdb                 ###                    DIL ILE  D-ISOLEUCINE
sed -i '/^ATOM/s/IIL/ILE/g' ${poi%.*}_r_prep1.pdb                 ###                    IIL ILE  ISO-ISOLEUCINE
sed -i '/^ATOM/s/ILE/ILE/g' ${poi%.*}_r_prep1.pdb                 ###                    ILE ILE
sed -i '/^ATOM/s/ILX/ILE/g' ${poi%.*}_r_prep1.pdb                 ###                    ILX ILE  45-DIHYDROXYISOLEUCINE
sed -i '/^ATOM/s/IML/ILE/g' ${poi%.*}_r_prep1.pdb                 ###                    IML ILE  N-METHYLATED
sed -i '/^ATOM/s/ALY/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
sed -i '/^ATOM/s/DLY/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
sed -i '/^ATOM/s/KCX/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
sed -i '/^ATOM/s/LLP/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
sed -i '/^ATOM/s/LLY/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
sed -i '/^ATOM/s/LYM/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
sed -i '/^ATOM/s/LYS/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                  LYS
sed -i '/^ATOM/s/LYZ/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
sed -i '/^ATOM/s/MLY/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
sed -i '/^ATOM/s/SHR/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
sed -i '/^ATOM/s/TRG/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
sed -i '/^ATOM/s/6CL/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    6CL LYS  6-CARBOXYLYSINE
sed -i '/^ATOM/s/ALY/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    ALY LYS  N(6)-ACETYLLYSINE
sed -i '/^ATOM/s/API/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    API LYS  26-DIAMINOPIMELIC ACID
sed -i '/^ATOM/s/APK/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    APK LYS
sed -i '/^ATOM/s/AZK/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    AZK LYS  (2S)-2-AMINO-6-TRIAZANYLHEXAN-1-OL
sed -i '/^ATOM/s/B3K/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    B3K LYS  (3S)-37-DIAMINOHEPTANOIC ACID
sed -i '/^ATOM/s/BLY/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    BLY LYS  LYSINE BORONIC ACID
sed -i '/^ATOM/s/C1X/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    C1X LYS  MODIFIED LYSINE
sed -i '/^ATOM/s/CLG/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CLG LYS
sed -i '/^ATOM/s/CLH/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CLH LYS
sed -i '/^ATOM/s/CYJ/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CYJ LYS  MODIFIED LYSINE
sed -i '/^ATOM/s/DLS/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    DLS LYS  DI-ACETYL-LYSINE
sed -i '/^ATOM/s/DLY/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    DLY LYS  D-LYSINE
sed -i '/^ATOM/s/DNL/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    DNL LYS  6-AMINO-HEXANAL
sed -i '/^ATOM/s/FHL/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    FHL LYS  MODIFIED LYSINE
sed -i '/^ATOM/s/GPL/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    GPL LYS  LYSINE GUANOSINE-5-MONOPHOSPHATE
sed -i '/^ATOM/s/IT1/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    IT1 LYS
sed -i '/^ATOM/s/KCX/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    KCX LYS  CARBAMOYLATED LYSINE
sed -i '/^ATOM/s/KGC/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    KGC LYS
sed -i '/^ATOM/s/KST/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    KST LYS  N~6~-(5-CARBOXY-3-THIENYL)-L-LYSINE
sed -i '/^ATOM/s/LA2/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LA2 LYS
sed -i '/^ATOM/s/LCK/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LCK LYS
sed -i '/^ATOM/s/LCX/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LCX LYS  CARBAMYLATED LYSINE
sed -i '/^ATOM/s/LDH/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LDH LYS  N~6~-ETHYL-L-LYSINE
sed -i '/^ATOM/s/LET/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LET LYS  ODIFIED LYSINE
sed -i '/^ATOM/s/LLP/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LLP LYS
sed -i '/^ATOM/s/LLY/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LLY LYS  NZ-(DICARBOXYMETHYL)LYSINE
sed -i '/^ATOM/s/LSO/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LSO LYS  MODIFIED LYSINE
sed -i '/^ATOM/s/LYM/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LYM LYS  DEOXY-METHYL-LYSINE
sed -i '/^ATOM/s/LYN/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LYN LYS  26-DIAMINO-HEXANOIC ACID AMIDE
sed -i '/^ATOM/s/LYP/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LYP LYS  N~6~-METHYL-N~6~-PROPYL-L-LYSINE
sed -i '/^ATOM/s/LYR/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LYR LYS  MODIFIED LYSINE
sed -i '/^ATOM/s/LYS/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LYS LYS
sed -i '/^ATOM/s/LYX/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LYX LYS  N-(2-COENZYME A)-PROPANOYL-LYSINE
sed -i '/^ATOM/s/LYZ/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LYZ LYS  5-HYDROXYLYSINE
sed -i '/^ATOM/s/M2L/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    M2L LYS
sed -i '/^ATOM/s/M3L/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    M3L LYS  N-TRIMETHYLLYSINE
sed -i '/^ATOM/s/MCL/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    MCL LYS  NZ-(1-CARBOXYETHYL)-LYSINE
sed -i '/^ATOM/s/MLY/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    MLY LYS  METHYLATED LYSINE
sed -i '/^ATOM/s/MLZ/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    MLZ LYS  N-METHYL-LYSINE
sed -i '/^ATOM/s/OBS/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    OBS LYS  MODIFIED LYSINE
sed -i '/^ATOM/s/SLZ/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    SLZ LYS  L-THIALYSINE
sed -i '/^ATOM/s/XX1/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    XX1 LYS  N~6~-7H-PURIN-6-YL-L-LYSINE
sed -i '/^ATOM/s/BUG/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
sed -i '/^ATOM/s/CLE/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
sed -i '/^ATOM/s/DLE/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
sed -i '/^ATOM/s/LEU/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                  LEU
sed -i '/^ATOM/s/MLE/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
sed -i '/^ATOM/s/NLE/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
sed -i '/^ATOM/s/NLN/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
sed -i '/^ATOM/s/NLP/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
sed -i '/^ATOM/s/1LU/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    1LU LEU  4-METHYL-PENTANOIC ACID-2-OXYL GROUP
sed -i '/^ATOM/s/2ML/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    2ML LEU  2-METHYLLEUCINE
sed -i '/^ATOM/s/BLE/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    BLE LEU  LEUCINE BORONIC ACID
sed -i '/^ATOM/s/BUG/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    BUG LEU  TERT-LEUCYL AMINE
sed -i '/^ATOM/s/CLE/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    CLE LEU  LEUCINE AMIDE
sed -i '/^ATOM/s/DCL/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    DCL LEU  2-AMINO-4-METHYL-PENTANYL GROUP
sed -i '/^ATOM/s/DLE/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    DLE LEU  D-LEUCINE
sed -i '/^ATOM/s/DNE/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    DNE LEU  D-NORLEUCINE
sed -i '/^ATOM/s/DNG/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    DNG LEU  N-FORMYL-D-NORLEUCINE
sed -i '/^ATOM/s/DNM/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    DNM LEU  D-N-METHYL NORLEUCINE
sed -i '/^ATOM/s/FLE/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    FLE LEU  FUROYL-LEUCINE
sed -i '/^ATOM/s/HLU/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    HLU LEU  BETA-HYDROXYLEUCINE
sed -i '/^ATOM/s/LED/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    LED LEU  POST-TRANSLATIONAL MODIFICATION
sed -i '/^ATOM/s/LEF/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    LEF LEU  2-5-FLUOROLEUCINE
sed -i '/^ATOM/s/LEU/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    LEU LEU
sed -i '/^ATOM/s/LNT/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    LNT LEU
sed -i '/^ATOM/s/MHL/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    MHL LEU  N-METHYLATED HYDROXY
sed -i '/^ATOM/s/MLE/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    MLE LEU  N-METHYLATED
sed -i '/^ATOM/s/MLL/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    MLL LEU  METHYL L-LEUCINATE
sed -i '/^ATOM/s/MNL/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    MNL LEU  4N-DIMETHYLNORLEUCINE
sed -i '/^ATOM/s/NLE/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    NLE LEU  NORLEUCINE
sed -i '/^ATOM/s/NLN/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    NLN LEU  NORLEUCINE AMIDE
sed -i '/^ATOM/s/NLO/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    NLO LEU  O-METHYL-L-NORLEUCINE
sed -i '/^ATOM/s/PLE/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    PLE LEU  LEUCINE PHOSPHINIC ACID
sed -i '/^ATOM/s/PPH/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    PPH LEU  PHENYLALANINE PHOSPHINIC ACID
sed -i '/^ATOM/s/CXM/MET/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS MET
sed -i '/^ATOM/s/FME/MET/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS MET
sed -i '/^ATOM/s/MET/MET/g' ${poi%.*}_r_prep1.pdb                 ###                  MET
sed -i '/^ATOM/s/MSE/MET/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS MET
sed -i '/^ATOM/s/OMT/MET/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS MET
sed -i '/^ATOM/s/AME/MET/g' ${poi%.*}_r_prep1.pdb                 ###                    AME MET  ACETYLATED METHIONINE
sed -i '/^ATOM/s/CXM/MET/g' ${poi%.*}_r_prep1.pdb                 ###                    CXM MET  N-CARBOXYMETHIONINE
sed -i '/^ATOM/s/ESC/MET/g' ${poi%.*}_r_prep1.pdb                 ###                    ESC MET  2-AMINO-4-ETHYL SULFANYL BUTYRIC ACID
sed -i '/^ATOM/s/FME/MET/g' ${poi%.*}_r_prep1.pdb                 ###                    FME MET  FORMYL-METHIONINE
sed -i '/^ATOM/s/FOR/MET/g' ${poi%.*}_r_prep1.pdb                 ###                    FOR MET
sed -i '/^ATOM/s/MET/MET/g' ${poi%.*}_r_prep1.pdb                 ###                    MET MET
sed -i '/^ATOM/s/MHO/MET/g' ${poi%.*}_r_prep1.pdb                 ###                    MHO MET  POST-TRANSLATIONAL MODIFICATION
sed -i '/^ATOM/s/MME/MET/g' ${poi%.*}_r_prep1.pdb                 ###                    MME MET  N-METHYL METHIONINE
sed -i '/^ATOM/s/MSE/MET/g' ${poi%.*}_r_prep1.pdb                 ###                    MSE MET  ELENOMETHIONINE
sed -i '/^ATOM/s/MSO/MET/g' ${poi%.*}_r_prep1.pdb                 ###                    MSO MET  METHIONINE SULFOXIDE
sed -i '/^ATOM/s/OMT/MET/g' ${poi%.*}_r_prep1.pdb                 ###                    OMT MET  METHIONINE SULFONE
sed -i '/^ATOM/s/SME/MET/g' ${poi%.*}_r_prep1.pdb                 ###                    SME MET  METHIONINE SULFOXIDE
sed -i '/^ATOM/s/ASN/ASN/g' ${poi%.*}_r_prep1.pdb                 ###                  ASN
sed -i '/^ATOM/s/MEN/ASN/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASN
sed -i '/^ATOM/s/AFA/ASN/g' ${poi%.*}_r_prep1.pdb                 ###                    AFA ASN  N-[7-METHYL-OCT-24-DIENOYL]ASPARAGINE
sed -i '/^ATOM/s/AHB/ASN/g' ${poi%.*}_r_prep1.pdb                 ###                    AHB ASN  BETA-HYDROXYASPARAGINE
sed -i '/^ATOM/s/ASN/ASN/g' ${poi%.*}_r_prep1.pdb                 ###                    ASN ASN
sed -i '/^ATOM/s/B3X/ASN/g' ${poi%.*}_r_prep1.pdb                 ###                    B3X ASN  (3S)-35-DIAMINO-5-OXOPENTANOIC ACID
sed -i '/^ATOM/s/DMH/ASN/g' ${poi%.*}_r_prep1.pdb                 ###                    DMH ASN  N4N4-DIMETHYL-ASPARAGINE
sed -i '/^ATOM/s/DSG/ASN/g' ${poi%.*}_r_prep1.pdb                 ###                    DSG ASN  D-ASPARAGINE
sed -i '/^ATOM/s/MEN/ASN/g' ${poi%.*}_r_prep1.pdb                 ###                    MEN ASN  GAMMA METHYL ASPARAGINE
sed -i '/^ATOM/s/DPR/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PRO
sed -i '/^ATOM/s/PRO/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                  PRO
sed -i '/^ATOM/s/1AB/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                    1AB PRO  14-DIDEOXY-14-IMINO-D-ARABINITOL
sed -i '/^ATOM/s/2MT/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                    2MT PRO
sed -i '/^ATOM/s/4FB/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                    4FB PRO  (4S)-4-FLUORO-L-PROLINE
sed -i '/^ATOM/s/DPL/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                    DPL PRO  4-OXOPROLINE
sed -i '/^ATOM/s/DPR/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                    DPR PRO  D-PROLINE
sed -i '/^ATOM/s/H5M/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                    H5M PRO  TRANS-3-HYDROXY-5-METHYLPROLINE
sed -i '/^ATOM/s/HY3/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                    HY3 PRO  3-HYDROXYPROLINE
sed -i '/^ATOM/s/HYP/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                    HYP PRO  4-HYDROXYPROLINE
sed -i '/^ATOM/s/LPD/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                    LPD PRO  L-PROLINAMIDE
sed -i '/^ATOM/s/P2Y/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                    P2Y PRO  (2S)-PYRROLIDIN-2-YLMETHYLAMINE
sed -i '/^ATOM/s/PCA/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                    PCA PRO  5-OXOPROLINE
sed -i '/^ATOM/s/POM/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                    POM PRO  CIS-5-METHYL-4-OXOPROLINE
sed -i '/^ATOM/s/PRO/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                    PRO PRO
sed -i '/^ATOM/s/PRS/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                    PRS PRO  THIOPROLINE
sed -i '/^ATOM/s/DGN/GLN/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLN
sed -i '/^ATOM/s/GLN/GLN/g' ${poi%.*}_r_prep1.pdb                 ###                  GLN
sed -i '/^ATOM/s/DGN/GLN/g' ${poi%.*}_r_prep1.pdb                 ###                    DGN GLN  D-GLUTAMINE
sed -i '/^ATOM/s/GHG/GLN/g' ${poi%.*}_r_prep1.pdb                 ###                    GHG GLN  GAMMA-HYDROXY-GLUTAMINE
sed -i '/^ATOM/s/GLH/GLN/g' ${poi%.*}_r_prep1.pdb                 ###                    GLH GLN
sed -i '/^ATOM/s/GLN/GLN/g' ${poi%.*}_r_prep1.pdb                 ###                    GLN GLN
sed -i '/^ATOM/s/MGN/GLN/g' ${poi%.*}_r_prep1.pdb                 ###                    MGN GLN  2-METHYL-GLUTAMINE
sed -i '/^ATOM/s/ACL/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
sed -i '/^ATOM/s/AGM/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
sed -i '/^ATOM/s/ARG/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                  ARG
sed -i '/^ATOM/s/ARM/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
sed -i '/^ATOM/s/DAR/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
sed -i '/^ATOM/s/HAR/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
sed -i '/^ATOM/s/HMR/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
sed -i '/^ATOM/s/2MR/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    2MR ARG  N3 N4-DIMETHYLARGININE
sed -i '/^ATOM/s/AAR/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    AAR ARG  ARGININEAMIDE
sed -i '/^ATOM/s/ACL/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    ACL ARG  DEOXY-CHLOROMETHYL-ARGININE
sed -i '/^ATOM/s/AGM/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    AGM ARG  4-METHYL-ARGININE
sed -i '/^ATOM/s/ALG/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    ALG ARG  GUANIDINOBUTYRYL GROUP
sed -i '/^ATOM/s/AR2/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    AR2 ARG  ARGINYL-BENZOTHIAZOLE-6-CARBOXYLIC ACID
sed -i '/^ATOM/s/ARG/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    ARG ARG
sed -i '/^ATOM/s/ARM/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    ARM ARG  DEOXY-METHYL-ARGININE
sed -i '/^ATOM/s/ARO/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    ARO ARG  C-GAMMA-HYDROXY ARGININE
sed -i '/^ATOM/s/BOR/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    BOR ARG
sed -i '/^ATOM/s/CIR/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    CIR ARG  CITRULLINE
sed -i '/^ATOM/s/DA2/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    DA2 ARG  MODIFIED ARGININE
sed -i '/^ATOM/s/DAR/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    DAR ARG  D-ARGININE
sed -i '/^ATOM/s/HMR/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    HMR ARG  BETA-HOMOARGININE
sed -i '/^ATOM/s/HRG/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    HRG ARG  L-HOMOARGININE
sed -i '/^ATOM/s/MAI/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    MAI ARG  DEOXO-METHYLARGININE
sed -i '/^ATOM/s/MGG/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    MGG ARG  MODIFIED D-ARGININE
sed -i '/^ATOM/s/NMM/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    NMM ARG  MODIFIED ARGININE
sed -i '/^ATOM/s/OPR/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    OPR ARG  C-(3-OXOPROPYL)ARGININE
sed -i '/^ATOM/s/ORQ/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    ORQ ARG  N~5~-ACETYL-L-ORNITHINE
sed -i '/^ATOM/s/TYZ/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    TYZ ARG  PARA ACETAMIDO BENZOIC ACID
sed -i '/^ATOM/s/DSN/SER/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
sed -i '/^ATOM/s/MIS/SER/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
sed -i '/^ATOM/s/OAS/SER/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
sed -i '/^ATOM/s/SAC/SER/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
sed -i '/^ATOM/s/SEL/SER/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
sed -i '/^ATOM/s/SEP/SER/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
sed -i '/^ATOM/s/SER/SER/g' ${poi%.*}_r_prep1.pdb                 ###                  SER
sed -i '/^ATOM/s/SET/SER/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
sed -i '/^ATOM/s/SVA/SER/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
sed -i '/^ATOM/s/B3S/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    B3S SER  (3R)-3-AMINO-4-HYDROXYBUTANOIC ACID
sed -i '/^ATOM/s/BG1/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    BG1 SER
sed -i '/^ATOM/s/DHL/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    DHL SER  POST-TRANSLATIONAL MODIFICATION
sed -i '/^ATOM/s/DSE/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    DSE SER  D-SERINE N-METHYLATED
sed -i '/^ATOM/s/DSN/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    DSN SER  D-SERINE
sed -i '/^ATOM/s/FGP/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    FGP SER
sed -i '/^ATOM/s/GVL/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    GVL SER  SERINE MODIFED WITH PHOSPHOPANTETHEINE
sed -i '/^ATOM/s/HSE/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    HSE SER  L-HOMOSERINE
sed -i '/^ATOM/s/HSL/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    HSL SER  HOMOSERINE LACTONE
sed -i '/^ATOM/s/MC1/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    MC1 SER  METHICILLIN ACYL-SERINE
sed -i '/^ATOM/s/MIS/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    MIS SER  MODIFIED SERINE
sed -i '/^ATOM/s/N10/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    N10 SER  O-[(HEXYLAMINO)CARBONYL]-L-SERINE
sed -i '/^ATOM/s/NC1/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    NC1 SER  NITROCEFIN ACYL-SERINE
sed -i '/^ATOM/s/OAS/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    OAS SER  O-ACETYLSERINE
sed -i '/^ATOM/s/OSE/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    OSE SER  O-SULFO-L-SERINE
sed -i '/^ATOM/s/PG1/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    PG1 SER  BENZYLPENICILLOYL-ACYLATED SERINE
sed -i '/^ATOM/s/PYR/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    PYR SER  CHEMICALLY MODIFIED
sed -i '/^ATOM/s/S1H/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    S1H SER  1-HEXADECANOSULFONYL-O-L-SERINE
sed -i '/^ATOM/s/SAC/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SAC SER  N-ACETYL-SERINE
sed -i '/^ATOM/s/SBD/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SBD SER
sed -i '/^ATOM/s/SBG/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SBG SER  MODIFIED SERINE
sed -i '/^ATOM/s/SBL/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SBL SER
sed -i '/^ATOM/s/SDP/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SDP SER
sed -i '/^ATOM/s/SEB/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SEB SER  O-BENZYLSULFONYL-SERINE
sed -i '/^ATOM/s/SEL/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SEL SER  2-AMINO-13-PROPANEDIOL
sed -i '/^ATOM/s/SEP/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SEP SER  E PHOSPHOSERINE
sed -i '/^ATOM/s/SER/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SER SER
sed -i '/^ATOM/s/SET/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SET SER  AMINOSERINE
sed -i '/^ATOM/s/SGB/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SGB SER  MODIFIED SERINE
sed -i '/^ATOM/s/SGR/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SGR SER  MODIFIED SERINE
sed -i '/^ATOM/s/SOY/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SOY SER  OXACILLOYL-ACYLATED SERINE
sed -i '/^ATOM/s/SUN/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SUN SER  TABUN CONJUGATED SERINE
sed -i '/^ATOM/s/SVA/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SVA SER  SERINE VANADATE
sed -i '/^ATOM/s/SVV/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SVV SER  MODIFIED SERINE
sed -i '/^ATOM/s/SVX/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SVX SER  MODIFIED SERINE
sed -i '/^ATOM/s/SVY/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SVY SER  MODIFIED SERINE
sed -i '/^ATOM/s/SVZ/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SVZ SER  MODIFIED SERINE
sed -i '/^ATOM/s/SXE/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SXE SER  MODIFIED SERINE
sed -i '/^ATOM/s/ALO/THR/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS THR
sed -i '/^ATOM/s/BMT/THR/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS THR
sed -i '/^ATOM/s/DTH/THR/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS THR
sed -i '/^ATOM/s/THR/THR/g' ${poi%.*}_r_prep1.pdb                 ###                  THR
sed -i '/^ATOM/s/TPO/THR/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS THR
sed -i '/^ATOM/s/AEI/THR/g' ${poi%.*}_r_prep1.pdb                 ###                    AEI THR  ACYLATED THR
sed -i '/^ATOM/s/ALO/THR/g' ${poi%.*}_r_prep1.pdb                 ###                    ALO THR  ALLO-THREONINE
sed -i '/^ATOM/s/BMT/THR/g' ${poi%.*}_r_prep1.pdb                 ###                    BMT THR
sed -i '/^ATOM/s/CRO/THR/g' ${poi%.*}_r_prep1.pdb                 ###                    CRO THR  CYCLIZED
sed -i '/^ATOM/s/CTH/THR/g' ${poi%.*}_r_prep1.pdb                 ###                    CTH THR  4-CHLOROTHREONINE
sed -i '/^ATOM/s/DTH/THR/g' ${poi%.*}_r_prep1.pdb                 ###                    DTH THR  D-THREONINE
sed -i '/^ATOM/s/OLT/THR/g' ${poi%.*}_r_prep1.pdb                 ###                    OLT THR  O-METHYL-L-THREONINE
sed -i '/^ATOM/s/TBM/THR/g' ${poi%.*}_r_prep1.pdb                 ###                    TBM THR
sed -i '/^ATOM/s/TH5/THR/g' ${poi%.*}_r_prep1.pdb                 ###                    TH5 THR  O-ACETYL-L-THREONINE
sed -i '/^ATOM/s/THC/THR/g' ${poi%.*}_r_prep1.pdb                 ###                    THC THR  N-METHYLCARBONYLTHREONINE
sed -i '/^ATOM/s/THR/THR/g' ${poi%.*}_r_prep1.pdb                 ###                    THR THR
sed -i '/^ATOM/s/TMD/THR/g' ${poi%.*}_r_prep1.pdb                 ###                    TMD THR  N-METHYLATED EPSILON C ALKYLATED
sed -i '/^ATOM/s/TPO/THR/g' ${poi%.*}_r_prep1.pdb                 ###                    TPO THR  HOSPHOTHREONINE
sed -i '/^ATOM/s/DIV/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS VAL
sed -i '/^ATOM/s/DVA/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS VAL
sed -i '/^ATOM/s/MVA/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS VAL
sed -i '/^ATOM/s/VAL/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                  VAL
sed -i '/^ATOM/s/B2V/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                    B2V VAL  VALINE BORONIC ACID
sed -i '/^ATOM/s/DIV/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                    DIV VAL  D-ISOVALINE
sed -i '/^ATOM/s/DVA/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                    DVA VAL  D-VALINE
sed -i '/^ATOM/s/MNV/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                    MNV VAL  N-METHYL-C-AMINO VALINE
sed -i '/^ATOM/s/MVA/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                    MVA VAL  N-METHYLATED
sed -i '/^ATOM/s/NVA/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                    NVA VAL  NORVALINE
sed -i '/^ATOM/s/VAD/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                    VAD VAL  DEAMINOHYDROXYVALINE
sed -i '/^ATOM/s/VAF/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                    VAF VAL  METHYLVALINE
sed -i '/^ATOM/s/VAL/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                    VAL VAL
sed -i '/^ATOM/s/VDL/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                    VDL VAL  (2R3R)-23-DIAMINOBUTANOIC ACID
sed -i '/^ATOM/s/VLL/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                    VLL VAL  (2S)-23-DIAMINOBUTANOIC ACID
sed -i '/^ATOM/s/VME/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                    VME VAL  O- METHYLVALINE
sed -i '/^ATOM/s/DTR/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
sed -i '/^ATOM/s/HTR/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
sed -i '/^ATOM/s/LTR/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
sed -i '/^ATOM/s/TPL/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
sed -i '/^ATOM/s/TRO/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
sed -i '/^ATOM/s/TRP/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                  TRP
sed -i '/^ATOM/s/BTR/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    BTR TRP  6-BROMO-TRYPTOPHAN
sed -i '/^ATOM/s/1TQ/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    1TQ TRP  6-(FORMYLAMINO)-7-HYDROXY-L-TRYPTOPHAN
sed -i '/^ATOM/s/23S/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    23S TRP  MODIFIED TRYPTOPHAN
sed -i '/^ATOM/s/32S/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    32S TRP  MODIFIED TRYPTOPHAN
sed -i '/^ATOM/s/32T/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    32T TRP  MODIFIED TRYPTOPHAN
sed -i '/^ATOM/s/4DP/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    4DP TRP
sed -i '/^ATOM/s/4FW/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    4FW TRP  4-FLUOROTRYPTOPHANE
sed -i '/^ATOM/s/4HT/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    4HT TRP  4-HYDROXYTRYPTOPHAN
sed -i '/^ATOM/s/4IN/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    4IN TRP  4-AMINO-L-TRYPTOPHAN
sed -i '/^ATOM/s/6CW/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    6CW TRP  6-CHLORO-L-TRYPTOPHAN
sed -i '/^ATOM/s/DTR/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    DTR TRP  D-TRYPTOPHAN
sed -i '/^ATOM/s/FTR/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    FTR TRP  FLUOROTRYPTOPHANE
sed -i '/^ATOM/s/HTR/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    HTR TRP  BETA-HYDROXYTRYPTOPHANE
sed -i '/^ATOM/s/PAT/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    PAT TRP  ALPHA-PHOSPHONO-TRYPTOPHAN
sed -i '/^ATOM/s/TOX/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    TOX TRP
sed -i '/^ATOM/s/TPL/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    TPL TRP  TRYTOPHANOL
sed -i '/^ATOM/s/TQQ/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    TQQ TRP
sed -i '/^ATOM/s/TRF/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    TRF TRP  N1-FORMYL-TRYPTOPHAN
sed -i '/^ATOM/s/TRN/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    TRN TRP  AZA-TRYPTOPHAN
sed -i '/^ATOM/s/TRO/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    TRO TRP  2-HYDROXY-TRYPTOPHAN
sed -i '/^ATOM/s/TRP/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    TRP TRP
sed -i '/^ATOM/s/TRQ/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    TRQ TRP
sed -i '/^ATOM/s/TRW/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    TRW TRP
sed -i '/^ATOM/s/TRX/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    TRX TRP  6-HYDROXYTRYPTOPHAN
sed -i '/^ATOM/s/TTQ/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    TTQ TRP  6-AMINO-7-HYDROXY-L-TRYPTOPHAN
sed -i '/^ATOM/s/DTY/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
sed -i '/^ATOM/s/IYR/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
sed -i '/^ATOM/s/PAQ/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
sed -i '/^ATOM/s/PTR/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
sed -i '/^ATOM/s/STY/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
sed -i '/^ATOM/s/TYB/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
sed -i '/^ATOM/s/TYQ/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
sed -i '/^ATOM/s/TYR/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
sed -i '/^ATOM/s/TYS/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                  TYR
sed -i '/^ATOM/s/TYY/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
sed -i '/^ATOM/s/1TY/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    1TY TYR
sed -i '/^ATOM/s/2TY/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    2TY TYR
sed -i '/^ATOM/s/3TY/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    3TY TYR  MODIFIED TYROSINE
sed -i '/^ATOM/s/B3Y/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    B3Y TYR
sed -i '/^ATOM/s/CRQ/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    CRQ TYR
sed -i '/^ATOM/s/DBY/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    DBY TYR  35 DIBROMOTYROSINE
sed -i '/^ATOM/s/DPQ/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    DPQ TYR  TYROSINE DERIVATIVE
sed -i '/^ATOM/s/DTY/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    DTY TYR  D-TYROSINE
sed -i '/^ATOM/s/ESB/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    ESB TYR
sed -i '/^ATOM/s/FLT/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    FLT TYR  FLUOROMALONYL TYROSINE
sed -i '/^ATOM/s/FTY/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    FTY TYR  DEOXY-DIFLUOROMETHELENE-PHOSPHOTYROSINE
sed -i '/^ATOM/s/IYR/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    IYR TYR  3-IODO-TYROSINE
sed -i '/^ATOM/s/MBQ/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    MBQ TYR
sed -i '/^ATOM/s/NIY/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    NIY TYR  META-NITRO-TYROSINE
sed -i '/^ATOM/s/NBQ/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    NBQ TYR
sed -i '/^ATOM/s/OTY/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    OTY TYR
sed -i '/^ATOM/s/PAQ/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    PAQ TYR  SEE REMARK 999
sed -i '/^ATOM/s/PTH/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    PTH TYR  METHYLENE-HYDROXY-PHOSPHOTYROSINE
sed -i '/^ATOM/s/PTM/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    PTM TYR  ALPHA-METHYL-O-PHOSPHOTYROSINE
sed -i '/^ATOM/s/PTR/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    PTR TYR  O-PHOSPHOTYROSINE
sed -i '/^ATOM/s/TCQ/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    TCQ TYR  MODIFIED TYROSINE
sed -i '/^ATOM/s/TTS/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    TTS TYR
sed -i '/^ATOM/s/TY2/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    TY2 TYR  3-AMINO-L-TYROSINE
sed -i '/^ATOM/s/TY3/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    TY3 TYR  3-HYDROXY-L-TYROSINE
sed -i '/^ATOM/s/TYB/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    TYB TYR  TYROSINAL
sed -i '/^ATOM/s/TYC/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    TYC TYR  L-TYROSINAMIDE
sed -i '/^ATOM/s/TYI/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    TYI TYR  35-DIIODOTYROSINE
sed -i '/^ATOM/s/TYN/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    TYN TYR  ADDUCT AT HYDROXY GROUP
sed -i '/^ATOM/s/TYO/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    TYO TYR
sed -i '/^ATOM/s/TYQ/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    TYQ TYR  AMINOQUINOL FORM OF TOPA QUINONONE
sed -i '/^ATOM/s/TYR/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    TYR TYR
sed -i '/^ATOM/s/TYS/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    TYS TYR  INE SULPHONATED TYROSINE
sed -i '/^ATOM/s/TYT/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    TYT TYR
sed -i '/^ATOM/s/TYY/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    TYY TYR  IMINOQUINONE FORM OF TOPA QUINONONE
sed -i '/^ATOM/s/YOF/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    YOF TYR  3-FLUOROTYROSINE


#################################################################################################################################################################

sed -i '/^ATOM/s/CSD/CYS/g' ${e3%.*}_r_prep1.pdb                 #3-SULFINOALANINE
sed -i '/^ATOM/s/HYP/PRO/g' ${e3%.*}_r_prep1.pdb                 #4-HYDROXYPROLINE
sed -i '/^ATOM/s/BMT/THR/g' ${e3%.*}_r_prep1.pdb                 #4-METHYL-4-[(E)-2-BUTENYL]-4,N-METHYL-THREONINE
sed -i '/^ATOM/s/5HP/GLU/g' ${e3%.*}_r_prep1.pdb                 #5-HYDROXYPROLINE
sed -i '/^ATOM/s/ABA/ALA/g' ${e3%.*}_r_prep1.pdb                 #ALPHA-AMINOBUTYRIC_ACID
sed -i '/^ATOM/s/AIB/ALA/g' ${e3%.*}_r_prep1.pdb                 #ALPHA-AMINOISOBUTYRIC_ACID
sed -i '/^ATOM/s/CSW/CYS/g' ${e3%.*}_r_prep1.pdb                 #CYSTEINE-S-DIOXIDE
sed -i '/^ATOM/s/OCS/CYS/g' ${e3%.*}_r_prep1.pdb                 #CYSTEINESULFONIC_ACID
sed -i '/^ATOM/s/DAL/ALA/g' ${e3%.*}_r_prep1.pdb                 #D-ALANINE
sed -i '/^ATOM/s/DAR/ARG/g' ${e3%.*}_r_prep1.pdb                 #D-ARGININE
sed -i '/^ATOM/s/DSG/ASN/g' ${e3%.*}_r_prep1.pdb                 #D-ASPARAGINE
sed -i '/^ATOM/s/DSP/ASP/g' ${e3%.*}_r_prep1.pdb                 #D-ASPARTATE
sed -i '/^ATOM/s/DCY/CYS/g' ${e3%.*}_r_prep1.pdb                 #D-CYSTEINE
sed -i '/^ATOM/s/CRO/CRO/g' ${e3%.*}_r_prep1.pdb                 #DECARBOXY(PARAHYDROXYBENZYLIDENE-IMIDAZOLIDINONE)THREONINE
sed -i '/^ATOM/s/DGL/GLU/g' ${e3%.*}_r_prep1.pdb                 #D-GLUTAMATE
sed -i '/^ATOM/s/DGN/GLN/g' ${e3%.*}_r_prep1.pdb                 #D-GLUTAMINE
sed -i '/^ATOM/s/DHI/HIS/g' ${e3%.*}_r_prep1.pdb                 #D-HISTIDINE
sed -i '/^ATOM/s/DIL/ILE/g' ${e3%.*}_r_prep1.pdb                 #D-ISOLEUCINE
sed -i '/^ATOM/s/DIV/VAL/g' ${e3%.*}_r_prep1.pdb                 #D-ISOVALINE
sed -i '/^ATOM/s/DLE/LEU/g' ${e3%.*}_r_prep1.pdb                 #D-LEUCINE
sed -i '/^ATOM/s/DLY/LYS/g' ${e3%.*}_r_prep1.pdb                 #D-LYSINE
sed -i '/^ATOM/s/DPN/PHE/g' ${e3%.*}_r_prep1.pdb                 #D-PHENYLALANINE
sed -i '/^ATOM/s/DPR/PRO/g' ${e3%.*}_r_prep1.pdb                 #D-PROLINE
sed -i '/^ATOM/s/DSN/SER/g' ${e3%.*}_r_prep1.pdb                 #D-SERINE
sed -i '/^ATOM/s/DTH/THR/g' ${e3%.*}_r_prep1.pdb                 #D-THREONINE
sed -i '/^ATOM/s/DTR/DTR/g' ${e3%.*}_r_prep1.pdb                 #D-TRYPTOPHANE
sed -i '/^ATOM/s/DTY/TYR/g' ${e3%.*}_r_prep1.pdb                 #D-TYROSINE
sed -i '/^ATOM/s/DVA/VAL/g' ${e3%.*}_r_prep1.pdb                 #D-VALINE
sed -i '/^ATOM/s/CGU/GLU/g' ${e3%.*}_r_prep1.pdb                 #GAMMA-CARBOXY-GLUTAMIC_ACID
sed -i '/^ATOM/s/KCX/LYS/g' ${e3%.*}_r_prep1.pdb                 #LYSINE_NZ-CARBOXYLIC_ACID
sed -i '/^ATOM/s/LLP/LYS/g' ${e3%.*}_r_prep1.pdb                 #LYSINE-PYRIDOXAL-5'-PHOSPHATE
sed -i '/^ATOM/s/CXM/MET/g' ${e3%.*}_r_prep1.pdb                 #N-CARBOXYMETHIONINE
sed -i '/^ATOM/s/FME/MET/g' ${e3%.*}_r_prep1.pdb                 #N-FORMYLMETHIONINE
sed -i '/^ATOM/s/MLE/LEU/g' ${e3%.*}_r_prep1.pdb                 #N-METHYLLEUCINE
sed -i '/^ATOM/s/MVA/VAL/g' ${e3%.*}_r_prep1.pdb                 #N-METHYLVALINE
sed -i '/^ATOM/s/NLE/LEU/g' ${e3%.*}_r_prep1.pdb                 #NORLEUCINE
sed -i '/^ATOM/s/PTR/TYR/g' ${e3%.*}_r_prep1.pdb                 #O-PHOSPHOTYROSINE
sed -i '/^ATOM/s/ORN/ALA/g' ${e3%.*}_r_prep1.pdb                 #ORNITHINE
sed -i '/^ATOM/s/SEP/SER/g' ${e3%.*}_r_prep1.pdb                 #PHOSPHOSERINE
sed -i '/^ATOM/s/TPO/THR/g' ${e3%.*}_r_prep1.pdb                 #PHOSPHOTHREONINE
sed -i '/^ATOM/s/PCA/GLU/g' ${e3%.*}_r_prep1.pdb                 #PYROGLUTAMIC_ACID
sed -i '/^ATOM/s/SAR/GLY/g' ${e3%.*}_r_prep1.pdb                 #SARCOSINE
sed -i '/^ATOM/s/CEA/CYS/g' ${e3%.*}_r_prep1.pdb                 #S-HYDROXY-CYSTEINE
sed -i '/^ATOM/s/CSO/CYS/g' ${e3%.*}_r_prep1.pdb                 #S-HYDROXYCYSTEINE
sed -i '/^ATOM/s/CSS/CYS/g' ${e3%.*}_r_prep1.pdb                 #S-MERCAPTOCYSTEINE
sed -i '/^ATOM/s/CSX/CYS/g' ${e3%.*}_r_prep1.pdb                 #S-OXY_CYSTEINE
sed -i '/^ATOM/s/CME/CYS/g' ${e3%.*}_r_prep1.pdb                 #S,S-(2-HYDROXYETHYL)THIOCYSTEINE
sed -i '/^ATOM/s/TYS/TYR/g' ${e3%.*}_r_prep1.pdb                 #SULFONATED_TYROSINE
sed -i '/^ATOM/s/TPQ/PHE/g' ${e3%.*}_r_prep1.pdb                 #TOPO-QUINONE
sed -i '/^ATOM/s/STY/TYR/g' ${e3%.*}_r_prep1.pdb                 #TYROSINE-O-SULPHONIC_ACID
sed -i '/^ATOM/s/CCS/CYS/g' ${e3%.*}_r_prep1.pdb                 #https://github.com/EvanBaugh/miscellaneous_scripts/blob/master/process_pdb.py
sed -i '/^ATOM/s/CALA/ALA/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CARG/ARG/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CASN/ASN/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CASP/ASP/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CCYS/CYS/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CCYX/CYX/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CGLN/GLN/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CGLU/GLU/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CGLY/GLY/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CHID/HID/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CHIE/HIE/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CHIP/HIP/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CHYP/HYP/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CILE/ILE/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CLEU/LEU/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CLYS/LYS/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CMET/MET/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CPHE/PHE/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CPRO/PRO/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CSER/SER/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CTHR/THR/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CTRP/TRP/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CTYR/TYR/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/CVAL/VAL/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NALA/ALA/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NARG/ARG/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NASN/ASN/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NASP/ASP/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NCYS/CYS/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NCYX/CYX/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NGLN/GLN/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NGLU/GLU/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NGLY/GLY/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NHID/HID/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NHIE/HIE/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NHIP/HIP/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NILE/ILE/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NLEU/LEU/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NLYS/LYS/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NMET/MET/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NPHE/PHE/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NPRO/PRO/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NSER/SER/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NTHR/THR/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NTRP/TRP/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NTYR/TYR/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/NVAL/VAL/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i '/^ATOM/s/DAS/ASP/g' ${e3%.*}_r_prep1.pdb                 #GUESS https://github.com/EvanBaugh/miscellaneous_scripts/blob/master/process_pdb.py
sed -i '/^ATOM/s/CAF/CYS/g' ${e3%.*}_r_prep1.pdb                 #GUESS https://github.com/EvanBaugh/miscellaneous_scripts/blob/master/process_pdb.py ## S-DIMETHYLARSINOYL-CYSTEINE
sed -i '/^ATOM/s/CAS/CYS/g' ${e3%.*}_r_prep1.pdb                 #GUESS https://github.com/EvanBaugh/miscellaneous_scripts/blob/master/process_pdb.py ## S-(DIMETHYLARSENIC)CYSTEINE
sed -i '/^ATOM/s/AIB/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/ALA/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  ALA
sed -i '/^ATOM/s/ALM/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/AYA/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/BNN/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/CHG/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/CSD/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/DAL/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/DHA/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/DNP/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/FLA/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/HAC/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/PRR/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/MAA/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/TIH/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/TPQ/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
sed -i '/^ATOM/s/0CS/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    0CS ALA  3-[(S)-HYDROPEROXYSULFINYL]-L-ALANINE
sed -i '/^ATOM/s/2BU/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    2BU ADE
sed -i '/^ATOM/s/2OP/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    2OP (2S  2-HYDROXYPROPANAL
sed -i '/^ATOM/s/4F3/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    4F3 ALA  CYCLIZED
sed -i '/^ATOM/s/AA4/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    AA4 ALA  2-AMINO-5-HYDROXYPENTANOIC ACID
sed -i '/^ATOM/s/ABA/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    ABA ALA  ALPHA-AMINOBUTYRIC ACID
sed -i '/^ATOM/s/AHO/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    AHO ALA  N-ACETYL-N-HYDROXY-L-ORNITHINE
sed -i '/^ATOM/s/AHP/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    AHP ALA  2-AMINO-HEPTANOIC ACID
sed -i '/^ATOM/s/AIB/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    AIB ALA  ALPHA-AMINOISOBUTYRIC ACID
sed -i '/^ATOM/s/ALA/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    ALA ALA
sed -i '/^ATOM/s/ALC/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    ALC ALA  2-AMINO-3-CYCLOHEXYL-PROPIONIC ACID
sed -i '/^ATOM/s/ALM/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    ALM ALA  1-METHYL-ALANINAL
sed -i '/^ATOM/s/ALN/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    ALN ALA  NAPHTHALEN-2-YL-3-ALANINE
sed -i '/^ATOM/s/ALS/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    ALS ALA  2-AMINO-3-OXO-4-SULFO-BUTYRIC ACID
sed -i '/^ATOM/s/ALT/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    ALT ALA  THIOALANINE
sed -i '/^ATOM/s/AP7/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    AP7 ADE
sed -i '/^ATOM/s/APH/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    APH ALA  P-AMIDINOPHENYL-3-ALANINE
sed -i '/^ATOM/s/AYA/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    AYA ALA  N-ACETYLALANINE
sed -i '/^ATOM/s/AYG/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    AYG ALA
sed -i '/^ATOM/s/B2A/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    B2A ALA  ALANINE BORONIC ACID
sed -i '/^ATOM/s/B3A/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    B3A ALA  (3S)-3-AMINOBUTANOIC ACID
sed -i '/^ATOM/s/BAL/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    BAL ALA  BETA-ALANINE
sed -i '/^ATOM/s/BNN/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    BNN ALA  ACETYL-P-AMIDINOPHENYLALANINE
sed -i '/^ATOM/s/C12/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    C12 ALA
sed -i '/^ATOM/s/C99/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    C99 ALA
sed -i '/^ATOM/s/CAB/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CAB ALA  4-CARBOXY-4-AMINOBUTANAL
sed -i '/^ATOM/s/CH6/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CH6 ALA
sed -i '/^ATOM/s/CH7/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CH7 ALA
sed -i '/^ATOM/s/CLB/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CLB ALA
sed -i '/^ATOM/s/CLD/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CLD ALA
sed -i '/^ATOM/s/CLV/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CLV ALA
sed -i '/^ATOM/s/CQR/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CQR ALA
sed -i '/^ATOM/s/CR2/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CR2 ALA  POST-TRANSLATIONAL MODIFICATION
sed -i '/^ATOM/s/CR5/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CR5 ALA
sed -i '/^ATOM/s/CR7/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CR7 ALA
sed -i '/^ATOM/s/CR8/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CR8 ALA
sed -i '/^ATOM/s/CRK/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CRK ALA
sed -i '/^ATOM/s/CRW/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CRW ALA
sed -i '/^ATOM/s/CRX/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CRX ALA
sed -i '/^ATOM/s/CSI/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CSI ALA
sed -i '/^ATOM/s/CSY/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CSY ALA  MODIFIED TYROSINE COMPLEX
sed -i '/^ATOM/s/CWR/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CWR ALA
sed -i '/^ATOM/s/DAB/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    DAB ALA  24-DIAMINOBUTYRIC ACID
sed -i '/^ATOM/s/DAL/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    DAL ALA  D-ALANINE
sed -i '/^ATOM/s/DAM/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    DAM ALA  N-METHYL-ALPHA-BETA-DEHYDROALANINE
sed -i '/^ATOM/s/DBU/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    DBU ALA  (2E)-2-AMINOBUT-2-ENOIC ACID
sed -i '/^ATOM/s/DBZ/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    DBZ ALA  3-(BENZOYLAMINO)-L-ALANINE
sed -i '/^ATOM/s/DHA/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    DHA ALA  2-AMINO-ACRYLIC ACID
sed -i '/^ATOM/s/DPP/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    DPP ALA  DIAMMINOPROPANOIC ACID
sed -i '/^ATOM/s/FGL/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    FGL ALA  2-AMINOPROPANEDIOIC ACID
sed -i '/^ATOM/s/HHK/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    HHK ALA  (2S)-28-DIAMINOOCTANOIC ACID
sed -i '/^ATOM/s/HMF/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    HMF ALA  2-AMINO-4-PHENYL-BUTYRIC ACID
sed -i '/^ATOM/s/IAM/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    IAM ALA  4-[(ISOPROPYLAMINO)METHYL]PHENYLALANINE
sed -i '/^ATOM/s/IGL/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    IGL ALA  ALPHA-AMINO-2-INDANACETIC ACID
sed -i '/^ATOM/s/KYN/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    KYN ALA  KYNURENINE
sed -i '/^ATOM/s/LAL/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    LAL ALA  NN-DIMETHYL-L-ALANINE
sed -i '/^ATOM/s/MAA/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    MAA ALA  N-METHYLALANINE
sed -i '/^ATOM/s/MDO/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    MDO ALA
sed -i '/^ATOM/s/MFC/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    MFC ALA  CYCLIZED
sed -i '/^ATOM/s/NAL/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    NAL ALA  BETA-(2-NAPHTHYL)-ALANINE
sed -i '/^ATOM/s/NAM/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    NAM ALA  NAM NAPTHYLAMINOALANINE
sed -i '/^ATOM/s/NCB/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    NCB ALA  CHEMICAL MODIFICATION
sed -i '/^ATOM/s/NRQ/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    NRQ ALA
sed -i '/^ATOM/s/NYC/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    NYC ALA
sed -i '/^ATOM/s/ORN/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    ORN ALA  ORNITHINE
sed -i '/^ATOM/s/PIA/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    PIA ALA  FUSION OF ALA 65 TYR 66 GLY 67
sed -i '/^ATOM/s/PRR/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    PRR ALA  3-(METHYL-PYRIDINIUM)ALANINE
sed -i '/^ATOM/s/PYA/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    PYA ALA  3-(110-PHENANTHROL-2-YL)-L-ALANINE
sed -i '/^ATOM/s/PYC/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    PYC ALA  PYRROLE-2-CARBOXYLATE
sed -i '/^ATOM/s/PYT/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    PYT ALA  MODIFIED ALANINE
sed -i '/^ATOM/s/RC7/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    RC7 ALA
sed -i '/^ATOM/s/SEC/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    SEC ALA  2-AMINO-3-SELENINO-PROPIONIC ACID
sed -i '/^ATOM/s/SIC/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    SIC ALA
sed -i '/^ATOM/s/SUI/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    SUI ALA
sed -i '/^ATOM/s/TIH/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    TIH ALA  BETA(2-THIENYL)ALANINE
sed -i '/^ATOM/s/TPQ/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    TPQ ALA  245-TRIHYDROXYPHENYLALANINE
sed -i '/^ATOM/s/UMA/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    UMA ALA
sed -i '/^ATOM/s/X9Q/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    X9Q ALA
sed -i '/^ATOM/s/XXY/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    XXY ALA
sed -i '/^ATOM/s/XYG/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    XYG ALA
sed -i '/^ATOM/s/BCS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/BUC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/C5C/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/C6C/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/CCS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/CEA/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/CME/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/CSO/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/CSP/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/CSS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/CSX/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/CSW/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/CY1/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/CY3/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/CYG/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/CYM/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/CYS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  CYS
sed -i '/^ATOM/s/CYQ/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/DCY/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/EFC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/OCS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/PEC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/PR3/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/SCH/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/SCS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/SCY/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/SHC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/SMC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/SOC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i '/^ATOM/s/5CS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    5CS CYS
sed -i '/^ATOM/s/AGT/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    AGT CYS  AGMATINE-CYSTEINE ADDUCT
sed -i '/^ATOM/s/BBC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    BBC CYS
sed -i '/^ATOM/s/BCS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    BCS CYS  BENZYLCYSTEINE
sed -i '/^ATOM/s/BCX/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    BCX CYS  BETA-3-CYSTEINE
sed -i '/^ATOM/s/BPE/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    BPE CYS
sed -i '/^ATOM/s/BUC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    BUC CYS  SS-BUTYLTHIOCYSTEINE
sed -i '/^ATOM/s/C3Y/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    C3Y CYS  MODIFIED CYSTEINE
sed -i '/^ATOM/s/C5C/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    C5C CYS  S-CYCLOPENTYL THIOCYSTEINE
sed -i '/^ATOM/s/C6C/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    C6C CYS  S-CYCLOHEXYL THIOCYSTEINE
sed -i '/^ATOM/s/CAF/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CAF CYS  S-DIMETHYLARSINOYL-CYSTEINE
sed -i '/^ATOM/s/CAS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CAS CYS  S-(DIMETHYLARSENIC)CYSTEINE
sed -i '/^ATOM/s/CCS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CCS CYS  CARBOXYMETHYLATED CYSTEINE
sed -i '/^ATOM/s/CME/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CME CYS  MODIFIED CYSTEINE
sed -i '/^ATOM/s/CML/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CML CYS
sed -i '/^ATOM/s/CMT/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CMT CYS  O-METHYLCYSTEINE
sed -i '/^ATOM/s/CS1/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CS1 CYS  S-(2-ANILINYL-SULFANYL)-CYSTEINE
sed -i '/^ATOM/s/CS3/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CS3 CYS
sed -i '/^ATOM/s/CS4/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CS4 CYS
sed -i '/^ATOM/s/CSA/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CSA CYS  S-ACETONYLCYSTEIN
sed -i '/^ATOM/s/CSB/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CSB CYS  CYS BOUND TO LEAD ION
sed -i '/^ATOM/s/CSD/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CSD CYS  3-SULFINOALANINE
sed -i '/^ATOM/s/CSE/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CSE CYS  SELENOCYSTEINE
sed -i '/^ATOM/s/CSO/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CSO CYS  INE S-HYDROXYCYSTEINE
sed -i '/^ATOM/s/CSR/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CSR CYS  S-ARSONOCYSTEINE
sed -i '/^ATOM/s/CSS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CSS CYS  13-THIAZOLE-4-CARBOXYLIC ACID
sed -i '/^ATOM/s/CSU/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CSU CYS  CYSTEINE-S-SULFONIC ACID
sed -i '/^ATOM/s/CSW/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CSW CYS  CYSTEINE-S-DIOXIDE
sed -i '/^ATOM/s/CSX/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CSX CYS  OXOCYSTEINE
sed -i '/^ATOM/s/CSZ/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CSZ CYS  S-SELANYL CYSTEINE
sed -i '/^ATOM/s/CY0/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CY0 CYS  MODIFIED CYSTEINE
sed -i '/^ATOM/s/CY1/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CY1 CYS  ACETAMIDOMETHYLCYSTEINE
sed -i '/^ATOM/s/CY3/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CY3 CYS  2-AMINO-3-MERCAPTO-PROPIONAMIDE
sed -i '/^ATOM/s/CY4/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CY4 CYS  S-BUTYRYL-CYSTEIN
sed -i '/^ATOM/s/CY7/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CY7 CYS  MODIFIED CYSTEINE
sed -i '/^ATOM/s/CYD/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CYD CYS
sed -i '/^ATOM/s/CYF/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CYF CYS  FLUORESCEIN LABELLED CYS380 (P14)
sed -i '/^ATOM/s/CYG/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CYG CYS
sed -i '/^ATOM/s/CYQ/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CYQ CYS
sed -i '/^ATOM/s/CYR/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CYR CYS
sed -i '/^ATOM/s/CYS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CYS CYS
sed -i '/^ATOM/s/CZ2/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CZ2 CYS  S-(DIHYDROXYARSINO)CYSTEINE
sed -i '/^ATOM/s/CZZ/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CZZ CYS  THIARSAHYDROXY-CYSTEINE
sed -i '/^ATOM/s/DCY/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    DCY CYS  D-CYSTEINE
sed -i '/^ATOM/s/DYS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    DYS CYS
sed -i '/^ATOM/s/EFC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    EFC CYS  SS-(2-FLUOROETHYL)THIOCYSTEINE
sed -i '/^ATOM/s/FOE/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    FOE CYS
sed -i '/^ATOM/s/GT9/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    GT9 CYS  SG ALKYLATED
sed -i '/^ATOM/s/GYC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    GYC CYS
sed -i '/^ATOM/s/HTI/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    HTI CYS
sed -i '/^ATOM/s/KOR/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    KOR CYS  MODIFIED CYSTEINE
sed -i '/^ATOM/s/M0H/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    M0H CYS  S-(HYDROXYMETHYL)-L-CYSTEINE
sed -i '/^ATOM/s/MCS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    MCS CYS  MALONYLCYSTEINE
sed -i '/^ATOM/s/NPH/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    NPH CYS
sed -i '/^ATOM/s/NYS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    NYS CYS
sed -i '/^ATOM/s/OCS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    OCS CYS  CYSTEINE SULFONIC ACID
sed -i '/^ATOM/s/OCY/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    OCY CYS  HYDROXYETHYLCYSTEINE
sed -i '/^ATOM/s/P1L/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    P1L CYS  S-PALMITOYL CYSTEINE
sed -i '/^ATOM/s/PBB/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    PBB CYS  S-(4-BROMOBENZYL)CYSTEINE
sed -i '/^ATOM/s/PEC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    PEC CYS  SS-PENTYLTHIOCYSTEINE
sed -i '/^ATOM/s/PR3/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    PR3 CYS  INE DTT-CYSTEINE
sed -i '/^ATOM/s/PYX/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    PYX CYS  S-[S-THIOPYRIDOXAMINYL]CYSTEINE
sed -i '/^ATOM/s/R1A/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    R1A CYS
sed -i '/^ATOM/s/R1B/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    R1B CYS
sed -i '/^ATOM/s/R1F/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    R1F CYS
sed -i '/^ATOM/s/R7A/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    R7A CYS
sed -i '/^ATOM/s/RCY/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    RCY CYS
sed -i '/^ATOM/s/SAH/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    SAH CYS  S-ADENOSYL-L-HOMOCYSTEINE
sed -i '/^ATOM/s/SC2/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    SC2 CYS  N-ACETYL-L-CYSTEINE
sed -i '/^ATOM/s/SCH/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    SCH CYS  S-METHYL THIOCYSTEINE GROUP
sed -i '/^ATOM/s/SCS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    SCS CYS  MODIFIED CYSTEINE
sed -i '/^ATOM/s/SCY/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    SCY CYS  CETYLATED CYSTEINE
sed -i '/^ATOM/s/SHC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    SHC CYS  S-HEXYLCYSTEINE
sed -i '/^ATOM/s/SMC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    SMC CYS  POST-TRANSLATIONAL MODIFICATION
sed -i '/^ATOM/s/SNC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    SNC CYS  S-NITROSO CYSTEINE
sed -i '/^ATOM/s/SOC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    SOC CYS  DIOXYSELENOCYSTEINE
sed -i '/^ATOM/s/TEE/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    TEE CYS  POST-TRANSLATIONAL MODIFICATION
sed -i '/^ATOM/s/TNB/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    TNB CYS  S-(236-TRINITROPHENYL)CYSTEINE
sed -i '/^ATOM/s/TYX/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    TYX CYS  S-(2-ANILINO-2-OXOETHYL)-L-CYSTEINE
sed -i '/^ATOM/s/YCM/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    YCM CYS  S-(2-AMINO-2-OXOETHYL)-L-CYSTEINE
sed -i '/^ATOM/s/2AS/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
sed -i '/^ATOM/s/ASA/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
sed -i '/^ATOM/s/ASB/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
sed -i '/^ATOM/s/ASK/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
sed -i '/^ATOM/s/ASL/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
sed -i '/^ATOM/s/ASP/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                  ASP
sed -i '/^ATOM/s/ASQ/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
sed -i '/^ATOM/s/BHD/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
sed -i '/^ATOM/s/DAS/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
sed -i '/^ATOM/s/DSP/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
sed -i '/^ATOM/s/3MD/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    3MD ASP  2S3S-3-METHYLASPARTIC ACID
sed -i '/^ATOM/s/A0A/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    A0A ASP  ASPARTYL-FORMYL MIXED ANHYDRIDE
sed -i '/^ATOM/s/ACB/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    ACB ASP  3-METHYL-ASPARTIC ACID
sed -i '/^ATOM/s/AKL/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    AKL ASP  3-AMINO-5-CHLORO-4-OXOPENTANOIC ACID
sed -i '/^ATOM/s/ASA/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    ASA ASP  ASPARTIC ALDEHYDE
sed -i '/^ATOM/s/ASB/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    ASB ASP  ASPARTIC ACID-4-CARBOXYETHYL ESTER
sed -i '/^ATOM/s/ASI/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    ASI ASP  L-ISO-ASPARTATE
sed -i '/^ATOM/s/ASK/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    ASK ASP  DEHYDROXYMETHYLASPARTIC ACID
sed -i '/^ATOM/s/ASL/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    ASL ASP  ASPARTIC ACID-4-CARBOXYETHYL ESTER
sed -i '/^ATOM/s/ASP/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    ASP ASP
sed -i '/^ATOM/s/B3D/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    B3D ASP  3-AMINOPENTANEDIOIC ACID
sed -i '/^ATOM/s/BFD/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    BFD ASP  ASPARTATE BERYLLIUM FLUORIDE
sed -i '/^ATOM/s/BHD/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    BHD ASP  BETA-HYDROXYASPARTIC ACID
sed -i '/^ATOM/s/DAS/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    DAS ASP  D-ASPARTIC ACID
sed -i '/^ATOM/s/DMK/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    DMK ASP  DIMETHYL ASPARTIC ACID
sed -i '/^ATOM/s/IAS/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    IAS ASP  ASPARTYL GROUP
sed -i '/^ATOM/s/OHS/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    OHS ASP  O-(CARBOXYSULFANYL)-4-OXO-L-HOMOSERINE
sed -i '/^ATOM/s/OXX/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    OXX ASP  OXALYL-ASPARTYL ANHYDRIDE
sed -i '/^ATOM/s/PHD/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    PHD ASP  2-AMINO-4-OXO-4-PHOSPHONOOXY-BUTYRIC ACID
sed -i '/^ATOM/s/SNN/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    SNN ASP  POST-TRANSLATIONAL MODIFICATION
sed -i '/^ATOM/s/5HP/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
sed -i '/^ATOM/s/CGU/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
sed -i '/^ATOM/s/DGL/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
sed -i '/^ATOM/s/GGL/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
sed -i '/^ATOM/s/GLU/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                  GLU
sed -i '/^ATOM/s/GMA/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
sed -i '/^ATOM/s/PCA/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
sed -i '/^ATOM/s/AB7/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                    AB7 GLU  ALPHA-AMINOBUTYRIC ACID
sed -i '/^ATOM/s/AR4/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                    AR4 GLU
sed -i '/^ATOM/s/B3E/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                    B3E GLU  (3S)-3-AMINOHEXANEDIOIC ACID
sed -i '/^ATOM/s/CGU/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                    CGU GLU  CARBOXYLATION OF THE CG ATOM
sed -i '/^ATOM/s/DGL/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                    DGL GLU  D-GLU
sed -i '/^ATOM/s/GLU/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                    GLU GLU
sed -i '/^ATOM/s/GMA/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                    GMA GLU  1-AMIDO-GLUTAMIC ACID
sed -i '/^ATOM/s/ILG/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                    ILG GLU  GLU LINKED TO NEXT RESIDUE VIA CG
sed -i '/^ATOM/s/LME/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                    LME GLU  (3R)-3-METHYL-L-GLUTAMIC ACID
sed -i '/^ATOM/s/MEG/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                    MEG GLU  (2S3R)-3-METHYL-GLUTAMIC ACID
sed -i '/^ATOM/s/DAH/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
sed -i '/^ATOM/s/DPN/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
sed -i '/^ATOM/s/HPQ/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
sed -i '/^ATOM/s/PHE/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                  PHE
sed -i '/^ATOM/s/PHI/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
sed -i '/^ATOM/s/PHL/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
sed -i '/^ATOM/s/1PA/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    1PA PHE  PHENYLMETHYLACETIC ACID ALANINE
sed -i '/^ATOM/s/23F/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    23F PHE  (2Z)-2-AMINO-3-PHENYLACRYLIC ACID
sed -i '/^ATOM/s/4PH/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    4PH PHE  4-METHYL-L-PHENYLALANINE
sed -i '/^ATOM/s/B2F/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    B2F PHE  PHENYLALANINE BORONIC ACID
sed -i '/^ATOM/s/BIF/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    BIF PHE
sed -i '/^ATOM/s/CHS/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    CHS PHE  4-AMINO-5-CYCLOHEXYL-3-HYDROXY-PENTANOIC AC
sed -i '/^ATOM/s/DAH/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    DAH PHE  34-DIHYDROXYDAHNYLALANINE
sed -i '/^ATOM/s/DPH/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    DPH PHE  DEAMINO-METHYL-PHENYLALANINE
sed -i '/^ATOM/s/DPN/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    DPN PHE  D-CONFIGURATION
sed -i '/^ATOM/s/FCL/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    FCL PHE  3-CHLORO-L-PHENYLALANINE
sed -i '/^ATOM/s/FOG/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    FOG PHE  PHENYLALANINOYL-[1-HYDROXY]-2-PROPYLENE
sed -i '/^ATOM/s/FRF/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    FRF PHE  PHE FOLLOWED BY REDUCED PHE
sed -i '/^ATOM/s/HPE/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    HPE PHE  HOMOPHENYLALANINE
sed -i '/^ATOM/s/HPH/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    HPH PHE  PHENYLALANINOL GROUP
sed -i '/^ATOM/s/HPQ/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    HPQ PHE  HOMOPHENYLALANINYLMETHANE
sed -i '/^ATOM/s/MEA/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    MEA PHE  N-METHYLPHENYLALANINE
sed -i '/^ATOM/s/MTY/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    MTY PHE  3-HYDROXYPHENYLALANINE
sed -i '/^ATOM/s/NFA/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    NFA PHE  MODIFIED PHENYLALANINE
sed -i '/^ATOM/s/PBF/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    PBF PHE  PARA-(BENZOYL)-PHENYLALANINE
sed -i '/^ATOM/s/PCS/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    PCS PHE  PHENYLALANYLMETHYLCHLORIDE
sed -i '/^ATOM/s/PF5/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    PF5 PHE  23456-PENTAFLUORO-L-PHENYLALANINE
sed -i '/^ATOM/s/PFF/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    PFF PHE  4-FLUORO-L-PHENYLALANINE
sed -i '/^ATOM/s/PHA/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    PHA PHE  PHENYLALANINAL
sed -i '/^ATOM/s/PHE/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    PHE PHE
sed -i '/^ATOM/s/PHI/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    PHI PHE  IODO-PHENYLALANINE
sed -i '/^ATOM/s/PHL/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    PHL PHE  L-PHENYLALANINOL
sed -i '/^ATOM/s/PHM/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    PHM PHE  PHENYLALANYLMETHANE
sed -i '/^ATOM/s/PM3/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    PM3 PHE
sed -i '/^ATOM/s/PPN/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    PPN PHE  THE LIGAND IS A PARA-NITRO-PHENYLALANINE
sed -i '/^ATOM/s/PRQ/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    PRQ PHE  PHENYLALANINE
sed -i '/^ATOM/s/PSA/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    PSA PHE
sed -i '/^ATOM/s/SMF/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    SMF PHE  4-SULFOMETHYL-L-PHENYLALANINE
sed -i '/^ATOM/s/GL3/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
sed -i '/^ATOM/s/GLY/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                  GLY
sed -i '/^ATOM/s/GLZ/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
sed -i '/^ATOM/s/GSC/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
sed -i '/^ATOM/s/MPQ/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
sed -i '/^ATOM/s/MSA/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
sed -i '/^ATOM/s/NMC/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
sed -i '/^ATOM/s/SAR/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
sed -i '/^ATOM/s/ACY/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    ACY GLY  POST-TRANSLATIONAL MODIFICATION
sed -i '/^ATOM/s/CHG/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    CHG GLY  CYCLOHEXYL GLYCINE
sed -i '/^ATOM/s/CHP/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    CHP GLY  3-CHLORO-4-HYDROXYPHENYLGLYCINE
sed -i '/^ATOM/s/GHP/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    GHP GLY  4-HYDROXYPHENYLGLYCINE
sed -i '/^ATOM/s/GL3/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    GL3 GLY  POST-TRANSLATIONAL MODIFICATION
sed -i '/^ATOM/s/GLY/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    GLY GLY
sed -i '/^ATOM/s/GLZ/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    GLZ GLY  AMINO-ACETALDEHYDE
sed -i '/^ATOM/s/GYS/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    GYS GLY
sed -i '/^ATOM/s/IPG/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    IPG GLY  N-ISOPROPYL GLYCINE
sed -i '/^ATOM/s/MEU/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    MEU GLY  O-METHYL-GLYCINE
sed -i '/^ATOM/s/MPQ/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    MPQ GLY  N-METHYL-ALPHA-PHENYL-GLYCINE
sed -i '/^ATOM/s/MSA/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    MSA GLY  (2-S-METHYL) SARCOSINE
sed -i '/^ATOM/s/NMC/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    NMC GLY  N-CYCLOPROPYLMETHYL GLYCINE
sed -i '/^ATOM/s/PG9/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    PG9 GLY  D-PHENYLGLYCINE
sed -i '/^ATOM/s/SAR/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    SAR GLY  SARCOSINE
sed -i '/^ATOM/s/SHP/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    SHP GLY  (4-HYDROXYMALTOSEPHENYL)GLYCINE
sed -i '/^ATOM/s/TBG/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    TBG GLY  T-BUTYL GLYCINE
sed -i '/^ATOM/s/3AH/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
sed -i '/^ATOM/s/DHI/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
sed -i '/^ATOM/s/HIC/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
sed -i '/^ATOM/s/HIS/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                  HIS
sed -i '/^ATOM/s/MHS/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
sed -i '/^ATOM/s/NEM/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
sed -i '/^ATOM/s/NEP/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
sed -i '/^ATOM/s/HID/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                  single delta N protonation
sed -i '/^ATOM/s/HIE/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                  single epsilon N protonation
sed -i '/^ATOM/s/3AH/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                    3AH HIS
sed -i '/^ATOM/s/DDE/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                    DDE HIS
sed -i '/^ATOM/s/DHI/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                    DHI HIS  D-HISTIDINE
sed -i '/^ATOM/s/HIA/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                    HIA HIS  L-HISTIDINE AMIDE
sed -i '/^ATOM/s/HIC/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                    HIC HIS  4-METHYL-HISTIDINE
sed -i '/^ATOM/s/HIP/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                    HIP HIS  ND1-PHOSPHONOHISTIDINE...or commonly used doubly protonated state
sed -i '/^ATOM/s/HIQ/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                    HIQ HIS  MODIFIED HISTIDINE
sed -i '/^ATOM/s/HIS/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                    HIS HIS
sed -i '/^ATOM/s/HSO/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                    HSO HIS  HISTIDINOL
sed -i '/^ATOM/s/MHS/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                    MHS HIS  1-N-METHYLHISTIDINE
sed -i '/^ATOM/s/NEP/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                    NEP HIS  N1-PHOSPHONOHISTIDINE
sed -i '/^ATOM/s/NZH/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                    NZH HIS
sed -i '/^ATOM/s/OHI/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                    OHI HIS  3-(2-OXO-2H-IMIDAZOL-4-YL)-L-ALANINE
sed -i '/^ATOM/s/PSH/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                    PSH HIS  1-THIOPHOSPHONO-L-HISTIDINE
sed -i '/^ATOM/s/DIL/ILE/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ILE
sed -i '/^ATOM/s/IIL/ILE/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ILE
sed -i '/^ATOM/s/ILE/ILE/g' ${e3%.*}_r_prep1.pdb                 ###                  ILE
sed -i '/^ATOM/s/B2I/ILE/g' ${e3%.*}_r_prep1.pdb                 ###                    B2I ILE  ISOLEUCINE BORONIC ACID
sed -i '/^ATOM/s/DIL/ILE/g' ${e3%.*}_r_prep1.pdb                 ###                    DIL ILE  D-ISOLEUCINE
sed -i '/^ATOM/s/IIL/ILE/g' ${e3%.*}_r_prep1.pdb                 ###                    IIL ILE  ISO-ISOLEUCINE
sed -i '/^ATOM/s/ILE/ILE/g' ${e3%.*}_r_prep1.pdb                 ###                    ILE ILE
sed -i '/^ATOM/s/ILX/ILE/g' ${e3%.*}_r_prep1.pdb                 ###                    ILX ILE  45-DIHYDROXYISOLEUCINE
sed -i '/^ATOM/s/IML/ILE/g' ${e3%.*}_r_prep1.pdb                 ###                    IML ILE  N-METHYLATED
sed -i '/^ATOM/s/ALY/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
sed -i '/^ATOM/s/DLY/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
sed -i '/^ATOM/s/KCX/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
sed -i '/^ATOM/s/LLP/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
sed -i '/^ATOM/s/LLY/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
sed -i '/^ATOM/s/LYM/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
sed -i '/^ATOM/s/LYS/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                  LYS
sed -i '/^ATOM/s/LYZ/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
sed -i '/^ATOM/s/MLY/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
sed -i '/^ATOM/s/SHR/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
sed -i '/^ATOM/s/TRG/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
sed -i '/^ATOM/s/6CL/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    6CL LYS  6-CARBOXYLYSINE
sed -i '/^ATOM/s/ALY/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    ALY LYS  N(6)-ACETYLLYSINE
sed -i '/^ATOM/s/API/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    API LYS  26-DIAMINOPIMELIC ACID
sed -i '/^ATOM/s/APK/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    APK LYS
sed -i '/^ATOM/s/AZK/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    AZK LYS  (2S)-2-AMINO-6-TRIAZANYLHEXAN-1-OL
sed -i '/^ATOM/s/B3K/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    B3K LYS  (3S)-37-DIAMINOHEPTANOIC ACID
sed -i '/^ATOM/s/BLY/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    BLY LYS  LYSINE BORONIC ACID
sed -i '/^ATOM/s/C1X/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    C1X LYS  MODIFIED LYSINE
sed -i '/^ATOM/s/CLG/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CLG LYS
sed -i '/^ATOM/s/CLH/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CLH LYS
sed -i '/^ATOM/s/CYJ/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CYJ LYS  MODIFIED LYSINE
sed -i '/^ATOM/s/DLS/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    DLS LYS  DI-ACETYL-LYSINE
sed -i '/^ATOM/s/DLY/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    DLY LYS  D-LYSINE
sed -i '/^ATOM/s/DNL/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    DNL LYS  6-AMINO-HEXANAL
sed -i '/^ATOM/s/FHL/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    FHL LYS  MODIFIED LYSINE
sed -i '/^ATOM/s/GPL/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    GPL LYS  LYSINE GUANOSINE-5-MONOPHOSPHATE
sed -i '/^ATOM/s/IT1/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    IT1 LYS
sed -i '/^ATOM/s/KCX/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    KCX LYS  CARBAMOYLATED LYSINE
sed -i '/^ATOM/s/KGC/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    KGC LYS
sed -i '/^ATOM/s/KST/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    KST LYS  N~6~-(5-CARBOXY-3-THIENYL)-L-LYSINE
sed -i '/^ATOM/s/LA2/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LA2 LYS
sed -i '/^ATOM/s/LCK/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LCK LYS
sed -i '/^ATOM/s/LCX/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LCX LYS  CARBAMYLATED LYSINE
sed -i '/^ATOM/s/LDH/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LDH LYS  N~6~-ETHYL-L-LYSINE
sed -i '/^ATOM/s/LET/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LET LYS  ODIFIED LYSINE
sed -i '/^ATOM/s/LLP/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LLP LYS
sed -i '/^ATOM/s/LLY/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LLY LYS  NZ-(DICARBOXYMETHYL)LYSINE
sed -i '/^ATOM/s/LSO/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LSO LYS  MODIFIED LYSINE
sed -i '/^ATOM/s/LYM/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LYM LYS  DEOXY-METHYL-LYSINE
sed -i '/^ATOM/s/LYN/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LYN LYS  26-DIAMINO-HEXANOIC ACID AMIDE
sed -i '/^ATOM/s/LYP/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LYP LYS  N~6~-METHYL-N~6~-PROPYL-L-LYSINE
sed -i '/^ATOM/s/LYR/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LYR LYS  MODIFIED LYSINE
sed -i '/^ATOM/s/LYS/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LYS LYS
sed -i '/^ATOM/s/LYX/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LYX LYS  N-(2-COENZYME A)-PROPANOYL-LYSINE
sed -i '/^ATOM/s/LYZ/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LYZ LYS  5-HYDROXYLYSINE
sed -i '/^ATOM/s/M2L/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    M2L LYS
sed -i '/^ATOM/s/M3L/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    M3L LYS  N-TRIMETHYLLYSINE
sed -i '/^ATOM/s/MCL/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    MCL LYS  NZ-(1-CARBOXYETHYL)-LYSINE
sed -i '/^ATOM/s/MLY/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    MLY LYS  METHYLATED LYSINE
sed -i '/^ATOM/s/MLZ/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    MLZ LYS  N-METHYL-LYSINE
sed -i '/^ATOM/s/OBS/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    OBS LYS  MODIFIED LYSINE
sed -i '/^ATOM/s/SLZ/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    SLZ LYS  L-THIALYSINE
sed -i '/^ATOM/s/XX1/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    XX1 LYS  N~6~-7H-PURIN-6-YL-L-LYSINE
sed -i '/^ATOM/s/BUG/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
sed -i '/^ATOM/s/CLE/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
sed -i '/^ATOM/s/DLE/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
sed -i '/^ATOM/s/LEU/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                  LEU
sed -i '/^ATOM/s/MLE/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
sed -i '/^ATOM/s/NLE/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
sed -i '/^ATOM/s/NLN/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
sed -i '/^ATOM/s/NLP/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
sed -i '/^ATOM/s/1LU/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    1LU LEU  4-METHYL-PENTANOIC ACID-2-OXYL GROUP
sed -i '/^ATOM/s/2ML/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    2ML LEU  2-METHYLLEUCINE
sed -i '/^ATOM/s/BLE/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    BLE LEU  LEUCINE BORONIC ACID
sed -i '/^ATOM/s/BUG/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    BUG LEU  TERT-LEUCYL AMINE
sed -i '/^ATOM/s/CLE/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    CLE LEU  LEUCINE AMIDE
sed -i '/^ATOM/s/DCL/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    DCL LEU  2-AMINO-4-METHYL-PENTANYL GROUP
sed -i '/^ATOM/s/DLE/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    DLE LEU  D-LEUCINE
sed -i '/^ATOM/s/DNE/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    DNE LEU  D-NORLEUCINE
sed -i '/^ATOM/s/DNG/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    DNG LEU  N-FORMYL-D-NORLEUCINE
sed -i '/^ATOM/s/DNM/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    DNM LEU  D-N-METHYL NORLEUCINE
sed -i '/^ATOM/s/FLE/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    FLE LEU  FUROYL-LEUCINE
sed -i '/^ATOM/s/HLU/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    HLU LEU  BETA-HYDROXYLEUCINE
sed -i '/^ATOM/s/LED/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    LED LEU  POST-TRANSLATIONAL MODIFICATION
sed -i '/^ATOM/s/LEF/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    LEF LEU  2-5-FLUOROLEUCINE
sed -i '/^ATOM/s/LEU/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    LEU LEU
sed -i '/^ATOM/s/LNT/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    LNT LEU
sed -i '/^ATOM/s/MHL/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    MHL LEU  N-METHYLATED HYDROXY
sed -i '/^ATOM/s/MLE/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    MLE LEU  N-METHYLATED
sed -i '/^ATOM/s/MLL/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    MLL LEU  METHYL L-LEUCINATE
sed -i '/^ATOM/s/MNL/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    MNL LEU  4N-DIMETHYLNORLEUCINE
sed -i '/^ATOM/s/NLE/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    NLE LEU  NORLEUCINE
sed -i '/^ATOM/s/NLN/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    NLN LEU  NORLEUCINE AMIDE
sed -i '/^ATOM/s/NLO/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    NLO LEU  O-METHYL-L-NORLEUCINE
sed -i '/^ATOM/s/PLE/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    PLE LEU  LEUCINE PHOSPHINIC ACID
sed -i '/^ATOM/s/PPH/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    PPH LEU  PHENYLALANINE PHOSPHINIC ACID
sed -i '/^ATOM/s/CXM/MET/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS MET
sed -i '/^ATOM/s/FME/MET/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS MET
sed -i '/^ATOM/s/MET/MET/g' ${e3%.*}_r_prep1.pdb                 ###                  MET
sed -i '/^ATOM/s/MSE/MET/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS MET
sed -i '/^ATOM/s/OMT/MET/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS MET
sed -i '/^ATOM/s/AME/MET/g' ${e3%.*}_r_prep1.pdb                 ###                    AME MET  ACETYLATED METHIONINE
sed -i '/^ATOM/s/CXM/MET/g' ${e3%.*}_r_prep1.pdb                 ###                    CXM MET  N-CARBOXYMETHIONINE
sed -i '/^ATOM/s/ESC/MET/g' ${e3%.*}_r_prep1.pdb                 ###                    ESC MET  2-AMINO-4-ETHYL SULFANYL BUTYRIC ACID
sed -i '/^ATOM/s/FME/MET/g' ${e3%.*}_r_prep1.pdb                 ###                    FME MET  FORMYL-METHIONINE
sed -i '/^ATOM/s/FOR/MET/g' ${e3%.*}_r_prep1.pdb                 ###                    FOR MET
sed -i '/^ATOM/s/MET/MET/g' ${e3%.*}_r_prep1.pdb                 ###                    MET MET
sed -i '/^ATOM/s/MHO/MET/g' ${e3%.*}_r_prep1.pdb                 ###                    MHO MET  POST-TRANSLATIONAL MODIFICATION
sed -i '/^ATOM/s/MME/MET/g' ${e3%.*}_r_prep1.pdb                 ###                    MME MET  N-METHYL METHIONINE
sed -i '/^ATOM/s/MSE/MET/g' ${e3%.*}_r_prep1.pdb                 ###                    MSE MET  ELENOMETHIONINE
sed -i '/^ATOM/s/MSO/MET/g' ${e3%.*}_r_prep1.pdb                 ###                    MSO MET  METHIONINE SULFOXIDE
sed -i '/^ATOM/s/OMT/MET/g' ${e3%.*}_r_prep1.pdb                 ###                    OMT MET  METHIONINE SULFONE
sed -i '/^ATOM/s/SME/MET/g' ${e3%.*}_r_prep1.pdb                 ###                    SME MET  METHIONINE SULFOXIDE
sed -i '/^ATOM/s/ASN/ASN/g' ${e3%.*}_r_prep1.pdb                 ###                  ASN
sed -i '/^ATOM/s/MEN/ASN/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASN
sed -i '/^ATOM/s/AFA/ASN/g' ${e3%.*}_r_prep1.pdb                 ###                    AFA ASN  N-[7-METHYL-OCT-24-DIENOYL]ASPARAGINE
sed -i '/^ATOM/s/AHB/ASN/g' ${e3%.*}_r_prep1.pdb                 ###                    AHB ASN  BETA-HYDROXYASPARAGINE
sed -i '/^ATOM/s/ASN/ASN/g' ${e3%.*}_r_prep1.pdb                 ###                    ASN ASN
sed -i '/^ATOM/s/B3X/ASN/g' ${e3%.*}_r_prep1.pdb                 ###                    B3X ASN  (3S)-35-DIAMINO-5-OXOPENTANOIC ACID
sed -i '/^ATOM/s/DMH/ASN/g' ${e3%.*}_r_prep1.pdb                 ###                    DMH ASN  N4N4-DIMETHYL-ASPARAGINE
sed -i '/^ATOM/s/DSG/ASN/g' ${e3%.*}_r_prep1.pdb                 ###                    DSG ASN  D-ASPARAGINE
sed -i '/^ATOM/s/MEN/ASN/g' ${e3%.*}_r_prep1.pdb                 ###                    MEN ASN  GAMMA METHYL ASPARAGINE
sed -i '/^ATOM/s/DPR/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PRO
sed -i '/^ATOM/s/PRO/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                  PRO
sed -i '/^ATOM/s/1AB/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                    1AB PRO  14-DIDEOXY-14-IMINO-D-ARABINITOL
sed -i '/^ATOM/s/2MT/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                    2MT PRO
sed -i '/^ATOM/s/4FB/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                    4FB PRO  (4S)-4-FLUORO-L-PROLINE
sed -i '/^ATOM/s/DPL/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                    DPL PRO  4-OXOPROLINE
sed -i '/^ATOM/s/DPR/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                    DPR PRO  D-PROLINE
sed -i '/^ATOM/s/H5M/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                    H5M PRO  TRANS-3-HYDROXY-5-METHYLPROLINE
sed -i '/^ATOM/s/HY3/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                    HY3 PRO  3-HYDROXYPROLINE
sed -i '/^ATOM/s/HYP/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                    HYP PRO  4-HYDROXYPROLINE
sed -i '/^ATOM/s/LPD/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                    LPD PRO  L-PROLINAMIDE
sed -i '/^ATOM/s/P2Y/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                    P2Y PRO  (2S)-PYRROLIDIN-2-YLMETHYLAMINE
sed -i '/^ATOM/s/PCA/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                    PCA PRO  5-OXOPROLINE
sed -i '/^ATOM/s/POM/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                    POM PRO  CIS-5-METHYL-4-OXOPROLINE
sed -i '/^ATOM/s/PRO/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                    PRO PRO
sed -i '/^ATOM/s/PRS/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                    PRS PRO  THIOPROLINE
sed -i '/^ATOM/s/DGN/GLN/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLN
sed -i '/^ATOM/s/GLN/GLN/g' ${e3%.*}_r_prep1.pdb                 ###                  GLN
sed -i '/^ATOM/s/DGN/GLN/g' ${e3%.*}_r_prep1.pdb                 ###                    DGN GLN  D-GLUTAMINE
sed -i '/^ATOM/s/GHG/GLN/g' ${e3%.*}_r_prep1.pdb                 ###                    GHG GLN  GAMMA-HYDROXY-GLUTAMINE
sed -i '/^ATOM/s/GLH/GLN/g' ${e3%.*}_r_prep1.pdb                 ###                    GLH GLN
sed -i '/^ATOM/s/GLN/GLN/g' ${e3%.*}_r_prep1.pdb                 ###                    GLN GLN
sed -i '/^ATOM/s/MGN/GLN/g' ${e3%.*}_r_prep1.pdb                 ###                    MGN GLN  2-METHYL-GLUTAMINE
sed -i '/^ATOM/s/ACL/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
sed -i '/^ATOM/s/AGM/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
sed -i '/^ATOM/s/ARG/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                  ARG
sed -i '/^ATOM/s/ARM/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
sed -i '/^ATOM/s/DAR/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
sed -i '/^ATOM/s/HAR/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
sed -i '/^ATOM/s/HMR/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
sed -i '/^ATOM/s/2MR/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    2MR ARG  N3 N4-DIMETHYLARGININE
sed -i '/^ATOM/s/AAR/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    AAR ARG  ARGININEAMIDE
sed -i '/^ATOM/s/ACL/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    ACL ARG  DEOXY-CHLOROMETHYL-ARGININE
sed -i '/^ATOM/s/AGM/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    AGM ARG  4-METHYL-ARGININE
sed -i '/^ATOM/s/ALG/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    ALG ARG  GUANIDINOBUTYRYL GROUP
sed -i '/^ATOM/s/AR2/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    AR2 ARG  ARGINYL-BENZOTHIAZOLE-6-CARBOXYLIC ACID
sed -i '/^ATOM/s/ARG/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    ARG ARG
sed -i '/^ATOM/s/ARM/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    ARM ARG  DEOXY-METHYL-ARGININE
sed -i '/^ATOM/s/ARO/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    ARO ARG  C-GAMMA-HYDROXY ARGININE
sed -i '/^ATOM/s/BOR/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    BOR ARG
sed -i '/^ATOM/s/CIR/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    CIR ARG  CITRULLINE
sed -i '/^ATOM/s/DA2/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    DA2 ARG  MODIFIED ARGININE
sed -i '/^ATOM/s/DAR/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    DAR ARG  D-ARGININE
sed -i '/^ATOM/s/HMR/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    HMR ARG  BETA-HOMOARGININE
sed -i '/^ATOM/s/HRG/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    HRG ARG  L-HOMOARGININE
sed -i '/^ATOM/s/MAI/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    MAI ARG  DEOXO-METHYLARGININE
sed -i '/^ATOM/s/MGG/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    MGG ARG  MODIFIED D-ARGININE
sed -i '/^ATOM/s/NMM/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    NMM ARG  MODIFIED ARGININE
sed -i '/^ATOM/s/OPR/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    OPR ARG  C-(3-OXOPROPYL)ARGININE
sed -i '/^ATOM/s/ORQ/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    ORQ ARG  N~5~-ACETYL-L-ORNITHINE
sed -i '/^ATOM/s/TYZ/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    TYZ ARG  PARA ACETAMIDO BENZOIC ACID
sed -i '/^ATOM/s/DSN/SER/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
sed -i '/^ATOM/s/MIS/SER/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
sed -i '/^ATOM/s/OAS/SER/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
sed -i '/^ATOM/s/SAC/SER/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
sed -i '/^ATOM/s/SEL/SER/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
sed -i '/^ATOM/s/SEP/SER/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
sed -i '/^ATOM/s/SER/SER/g' ${e3%.*}_r_prep1.pdb                 ###                  SER
sed -i '/^ATOM/s/SET/SER/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
sed -i '/^ATOM/s/SVA/SER/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
sed -i '/^ATOM/s/B3S/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    B3S SER  (3R)-3-AMINO-4-HYDROXYBUTANOIC ACID
sed -i '/^ATOM/s/BG1/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    BG1 SER
sed -i '/^ATOM/s/DHL/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    DHL SER  POST-TRANSLATIONAL MODIFICATION
sed -i '/^ATOM/s/DSE/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    DSE SER  D-SERINE N-METHYLATED
sed -i '/^ATOM/s/DSN/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    DSN SER  D-SERINE
sed -i '/^ATOM/s/FGP/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    FGP SER
sed -i '/^ATOM/s/GVL/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    GVL SER  SERINE MODIFED WITH PHOSPHOPANTETHEINE
sed -i '/^ATOM/s/HSE/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    HSE SER  L-HOMOSERINE
sed -i '/^ATOM/s/HSL/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    HSL SER  HOMOSERINE LACTONE
sed -i '/^ATOM/s/MC1/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    MC1 SER  METHICILLIN ACYL-SERINE
sed -i '/^ATOM/s/MIS/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    MIS SER  MODIFIED SERINE
sed -i '/^ATOM/s/N10/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    N10 SER  O-[(HEXYLAMINO)CARBONYL]-L-SERINE
sed -i '/^ATOM/s/NC1/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    NC1 SER  NITROCEFIN ACYL-SERINE
sed -i '/^ATOM/s/OAS/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    OAS SER  O-ACETYLSERINE
sed -i '/^ATOM/s/OSE/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    OSE SER  O-SULFO-L-SERINE
sed -i '/^ATOM/s/PG1/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    PG1 SER  BENZYLPENICILLOYL-ACYLATED SERINE
sed -i '/^ATOM/s/PYR/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    PYR SER  CHEMICALLY MODIFIED
sed -i '/^ATOM/s/S1H/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    S1H SER  1-HEXADECANOSULFONYL-O-L-SERINE
sed -i '/^ATOM/s/SAC/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SAC SER  N-ACETYL-SERINE
sed -i '/^ATOM/s/SBD/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SBD SER
sed -i '/^ATOM/s/SBG/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SBG SER  MODIFIED SERINE
sed -i '/^ATOM/s/SBL/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SBL SER
sed -i '/^ATOM/s/SDP/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SDP SER
sed -i '/^ATOM/s/SEB/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SEB SER  O-BENZYLSULFONYL-SERINE
sed -i '/^ATOM/s/SEL/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SEL SER  2-AMINO-13-PROPANEDIOL
sed -i '/^ATOM/s/SEP/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SEP SER  E PHOSPHOSERINE
sed -i '/^ATOM/s/SER/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SER SER
sed -i '/^ATOM/s/SET/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SET SER  AMINOSERINE
sed -i '/^ATOM/s/SGB/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SGB SER  MODIFIED SERINE
sed -i '/^ATOM/s/SGR/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SGR SER  MODIFIED SERINE
sed -i '/^ATOM/s/SOY/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SOY SER  OXACILLOYL-ACYLATED SERINE
sed -i '/^ATOM/s/SUN/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SUN SER  TABUN CONJUGATED SERINE
sed -i '/^ATOM/s/SVA/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SVA SER  SERINE VANADATE
sed -i '/^ATOM/s/SVV/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SVV SER  MODIFIED SERINE
sed -i '/^ATOM/s/SVX/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SVX SER  MODIFIED SERINE
sed -i '/^ATOM/s/SVY/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SVY SER  MODIFIED SERINE
sed -i '/^ATOM/s/SVZ/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SVZ SER  MODIFIED SERINE
sed -i '/^ATOM/s/SXE/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SXE SER  MODIFIED SERINE
sed -i '/^ATOM/s/ALO/THR/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS THR
sed -i '/^ATOM/s/BMT/THR/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS THR
sed -i '/^ATOM/s/DTH/THR/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS THR
sed -i '/^ATOM/s/THR/THR/g' ${e3%.*}_r_prep1.pdb                 ###                  THR
sed -i '/^ATOM/s/TPO/THR/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS THR
sed -i '/^ATOM/s/AEI/THR/g' ${e3%.*}_r_prep1.pdb                 ###                    AEI THR  ACYLATED THR
sed -i '/^ATOM/s/ALO/THR/g' ${e3%.*}_r_prep1.pdb                 ###                    ALO THR  ALLO-THREONINE
sed -i '/^ATOM/s/BMT/THR/g' ${e3%.*}_r_prep1.pdb                 ###                    BMT THR
sed -i '/^ATOM/s/CRO/THR/g' ${e3%.*}_r_prep1.pdb                 ###                    CRO THR  CYCLIZED
sed -i '/^ATOM/s/CTH/THR/g' ${e3%.*}_r_prep1.pdb                 ###                    CTH THR  4-CHLOROTHREONINE
sed -i '/^ATOM/s/DTH/THR/g' ${e3%.*}_r_prep1.pdb                 ###                    DTH THR  D-THREONINE
sed -i '/^ATOM/s/OLT/THR/g' ${e3%.*}_r_prep1.pdb                 ###                    OLT THR  O-METHYL-L-THREONINE
sed -i '/^ATOM/s/TBM/THR/g' ${e3%.*}_r_prep1.pdb                 ###                    TBM THR
sed -i '/^ATOM/s/TH5/THR/g' ${e3%.*}_r_prep1.pdb                 ###                    TH5 THR  O-ACETYL-L-THREONINE
sed -i '/^ATOM/s/THC/THR/g' ${e3%.*}_r_prep1.pdb                 ###                    THC THR  N-METHYLCARBONYLTHREONINE
sed -i '/^ATOM/s/THR/THR/g' ${e3%.*}_r_prep1.pdb                 ###                    THR THR
sed -i '/^ATOM/s/TMD/THR/g' ${e3%.*}_r_prep1.pdb                 ###                    TMD THR  N-METHYLATED EPSILON C ALKYLATED
sed -i '/^ATOM/s/TPO/THR/g' ${e3%.*}_r_prep1.pdb                 ###                    TPO THR  HOSPHOTHREONINE
sed -i '/^ATOM/s/DIV/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS VAL
sed -i '/^ATOM/s/DVA/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS VAL
sed -i '/^ATOM/s/MVA/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS VAL
sed -i '/^ATOM/s/VAL/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                  VAL
sed -i '/^ATOM/s/B2V/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                    B2V VAL  VALINE BORONIC ACID
sed -i '/^ATOM/s/DIV/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                    DIV VAL  D-ISOVALINE
sed -i '/^ATOM/s/DVA/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                    DVA VAL  D-VALINE
sed -i '/^ATOM/s/MNV/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                    MNV VAL  N-METHYL-C-AMINO VALINE
sed -i '/^ATOM/s/MVA/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                    MVA VAL  N-METHYLATED
sed -i '/^ATOM/s/NVA/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                    NVA VAL  NORVALINE
sed -i '/^ATOM/s/VAD/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                    VAD VAL  DEAMINOHYDROXYVALINE
sed -i '/^ATOM/s/VAF/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                    VAF VAL  METHYLVALINE
sed -i '/^ATOM/s/VAL/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                    VAL VAL
sed -i '/^ATOM/s/VDL/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                    VDL VAL  (2R3R)-23-DIAMINOBUTANOIC ACID
sed -i '/^ATOM/s/VLL/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                    VLL VAL  (2S)-23-DIAMINOBUTANOIC ACID
sed -i '/^ATOM/s/VME/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                    VME VAL  O- METHYLVALINE
sed -i '/^ATOM/s/DTR/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
sed -i '/^ATOM/s/HTR/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
sed -i '/^ATOM/s/LTR/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
sed -i '/^ATOM/s/TPL/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
sed -i '/^ATOM/s/TRO/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
sed -i '/^ATOM/s/TRP/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                  TRP
sed -i '/^ATOM/s/BTR/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    BTR TRP  6-BROMO-TRYPTOPHAN
sed -i '/^ATOM/s/1TQ/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    1TQ TRP  6-(FORMYLAMINO)-7-HYDROXY-L-TRYPTOPHAN
sed -i '/^ATOM/s/23S/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    23S TRP  MODIFIED TRYPTOPHAN
sed -i '/^ATOM/s/32S/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    32S TRP  MODIFIED TRYPTOPHAN
sed -i '/^ATOM/s/32T/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    32T TRP  MODIFIED TRYPTOPHAN
sed -i '/^ATOM/s/4DP/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    4DP TRP
sed -i '/^ATOM/s/4FW/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    4FW TRP  4-FLUOROTRYPTOPHANE
sed -i '/^ATOM/s/4HT/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    4HT TRP  4-HYDROXYTRYPTOPHAN
sed -i '/^ATOM/s/4IN/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    4IN TRP  4-AMINO-L-TRYPTOPHAN
sed -i '/^ATOM/s/6CW/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    6CW TRP  6-CHLORO-L-TRYPTOPHAN
sed -i '/^ATOM/s/DTR/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    DTR TRP  D-TRYPTOPHAN
sed -i '/^ATOM/s/FTR/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    FTR TRP  FLUOROTRYPTOPHANE
sed -i '/^ATOM/s/HTR/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    HTR TRP  BETA-HYDROXYTRYPTOPHANE
sed -i '/^ATOM/s/PAT/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    PAT TRP  ALPHA-PHOSPHONO-TRYPTOPHAN
sed -i '/^ATOM/s/TOX/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    TOX TRP
sed -i '/^ATOM/s/TPL/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    TPL TRP  TRYTOPHANOL
sed -i '/^ATOM/s/TQQ/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    TQQ TRP
sed -i '/^ATOM/s/TRF/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    TRF TRP  N1-FORMYL-TRYPTOPHAN
sed -i '/^ATOM/s/TRN/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    TRN TRP  AZA-TRYPTOPHAN
sed -i '/^ATOM/s/TRO/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    TRO TRP  2-HYDROXY-TRYPTOPHAN
sed -i '/^ATOM/s/TRP/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    TRP TRP
sed -i '/^ATOM/s/TRQ/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    TRQ TRP
sed -i '/^ATOM/s/TRW/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    TRW TRP
sed -i '/^ATOM/s/TRX/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    TRX TRP  6-HYDROXYTRYPTOPHAN
sed -i '/^ATOM/s/TTQ/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    TTQ TRP  6-AMINO-7-HYDROXY-L-TRYPTOPHAN
sed -i '/^ATOM/s/DTY/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
sed -i '/^ATOM/s/IYR/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
sed -i '/^ATOM/s/PAQ/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
sed -i '/^ATOM/s/PTR/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
sed -i '/^ATOM/s/STY/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
sed -i '/^ATOM/s/TYB/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
sed -i '/^ATOM/s/TYQ/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
sed -i '/^ATOM/s/TYR/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
sed -i '/^ATOM/s/TYS/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                  TYR
sed -i '/^ATOM/s/TYY/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
sed -i '/^ATOM/s/1TY/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    1TY TYR
sed -i '/^ATOM/s/2TY/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    2TY TYR
sed -i '/^ATOM/s/3TY/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    3TY TYR  MODIFIED TYROSINE
sed -i '/^ATOM/s/B3Y/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    B3Y TYR
sed -i '/^ATOM/s/CRQ/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    CRQ TYR
sed -i '/^ATOM/s/DBY/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    DBY TYR  35 DIBROMOTYROSINE
sed -i '/^ATOM/s/DPQ/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    DPQ TYR  TYROSINE DERIVATIVE
sed -i '/^ATOM/s/DTY/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    DTY TYR  D-TYROSINE
sed -i '/^ATOM/s/ESB/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    ESB TYR
sed -i '/^ATOM/s/FLT/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    FLT TYR  FLUOROMALONYL TYROSINE
sed -i '/^ATOM/s/FTY/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    FTY TYR  DEOXY-DIFLUOROMETHELENE-PHOSPHOTYROSINE
sed -i '/^ATOM/s/IYR/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    IYR TYR  3-IODO-TYROSINE
sed -i '/^ATOM/s/MBQ/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    MBQ TYR
sed -i '/^ATOM/s/NIY/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    NIY TYR  META-NITRO-TYROSINE
sed -i '/^ATOM/s/NBQ/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    NBQ TYR
sed -i '/^ATOM/s/OTY/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    OTY TYR
sed -i '/^ATOM/s/PAQ/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    PAQ TYR  SEE REMARK 999
sed -i '/^ATOM/s/PTH/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    PTH TYR  METHYLENE-HYDROXY-PHOSPHOTYROSINE
sed -i '/^ATOM/s/PTM/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    PTM TYR  ALPHA-METHYL-O-PHOSPHOTYROSINE
sed -i '/^ATOM/s/PTR/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    PTR TYR  O-PHOSPHOTYROSINE
sed -i '/^ATOM/s/TCQ/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    TCQ TYR  MODIFIED TYROSINE
sed -i '/^ATOM/s/TTS/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    TTS TYR
sed -i '/^ATOM/s/TY2/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    TY2 TYR  3-AMINO-L-TYROSINE
sed -i '/^ATOM/s/TY3/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    TY3 TYR  3-HYDROXY-L-TYROSINE
sed -i '/^ATOM/s/TYB/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    TYB TYR  TYROSINAL
sed -i '/^ATOM/s/TYC/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    TYC TYR  L-TYROSINAMIDE
sed -i '/^ATOM/s/TYI/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    TYI TYR  35-DIIODOTYROSINE
sed -i '/^ATOM/s/TYN/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    TYN TYR  ADDUCT AT HYDROXY GROUP
sed -i '/^ATOM/s/TYO/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    TYO TYR
sed -i '/^ATOM/s/TYQ/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    TYQ TYR  AMINOQUINOL FORM OF TOPA QUINONONE
sed -i '/^ATOM/s/TYR/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    TYR TYR
sed -i '/^ATOM/s/TYS/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    TYS TYR  INE SULPHONATED TYROSINE
sed -i '/^ATOM/s/TYT/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    TYT TYR
sed -i '/^ATOM/s/TYY/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    TYY TYR  IMINOQUINONE FORM OF TOPA QUINONONE
sed -i '/^ATOM/s/YOF/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    YOF TYR  3-FLUOROTYROSINE

#############################################################################################################################################################


pdb4amber -i ${poi%.*}_r_prep1.pdb -o ${poi%.*}_r_prep11.pdb
pdb4amber -i ${e3%.*}_r_prep1.pdb -o ${e3%.*}_r_prep11.pdb



sed -i '/UNL/!d' ${poi%.*}_ligand.pdb

sed -i '/UNL/!d' ${e3%.*}_ligand.pdb


sed -i '/HETATM.*ZN\|ZN.*HETATM/d' ${poi%.*}_ligand.pdb
sed -i '/HETATM.*ZN\|ZN.*HETATM/d' ${e3%.*}_ligand.pdb
sed -i '/HETATM.*NI\|NI.*HETATM/d' ${poi%.*}_ligand.pdb
sed -i '/HETATM.*NI\|NI.*HETATM/d' ${e3%.*}_ligand.pdb
sed -i '/HETATM.*MG\|MG.*HETATM/d' ${poi%.*}_ligand.pdb
sed -i '/HETATM.*MG\|MG.*HETATM/d' ${e3%.*}_ligand.pdb
sed -i '/HETATM.*CA\|CA.*HETATM/d' ${poi%.*}_ligand.pdb
sed -i '/HETATM.*CA\|CA.*HETATM/d' ${e3%.*}_ligand.pdb

sed -i 's/ATOM  /HETATM/g' ${poi%.*}_ligand.pdb
sed -i 's/ATOM  /HETATM/g' ${e3%.*}_ligand.pdb
sed -i 's/CL/Cl/g' ${poi%.*}_ligand.pdb
sed -i 's/CL/Cl/g' ${e3%.*}_ligand.pdb
sed -i 's/BR/Br/g' ${poi%.*}_ligand.pdb
sed -i 's/BR/Br/g' ${e3%.*}_ligand.pdb


#modified out oprot
#pdb4amber -i ${poi%.*}_ligand.pdb -o ${poi%.*}_ligand_oprot_prep1.pdb --reduce
#pdb4amber -i ${e3%.*}_ligand.pdb -o ${e3%.*}_ligand_oprot_prep1.pdb --reduce

#### EDITED
antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo mol2 -o ${poi%.*}_ligand_oprot_prep1_gas.mol2 -c gas -pf y

antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -dr n
#antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc 0 -pf y
else
	:
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -1 -pf y
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +1 -pf y
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -2 -pf y
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +2 -pf y
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -3 -pf y
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +3 -pf y
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +4 -pf y
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -4 -pf y
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -5 -pf y
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +5 -pf y
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
	#pdb4amber -i ${poi%.*}_ligand.pdb -o ${poi%.*}_ligand_oprot_prep1.pdb --reduce
	#antechamber -fi pdb -i ${poi%.*}_ligand_oprot_prep1.pdb -fo mol2 -o ${poi%.*}_ligand_oprot_prep1_gas.mol2 -c gas -pf y
	antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -dr n
else
	:
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc 0 -pf y -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -1 -pf y -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +1 -pf y -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -2 -pf y -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +2 -pf y -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -3 -pf y -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +3 -pf y -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +4 -pf y -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -4 -pf y -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -5 -pf y -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +5 -pf y -dr n
else
        :
fi


if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas.mol2 ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas.mol2 ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 0
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 1
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc -1
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 2
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc -2
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 3
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc -3
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas.mol2 ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas.mol2 ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 0 -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 1 -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc -1 -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 2 -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc -2 -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 3 -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc -3 -dr n
else
        :
fi

antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo mol2 -o ${e3%.*}_ligand_oprot_prep1_gas.mol2 -c gas -pf y
#antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -dr n
antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc 0 -pf y
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -1 -pf y
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +1 -pf y
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -2 -pf y
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +2 -pf y
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -3 -pf y
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +3 -pf y
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +4 -pf y
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -4 -pf y
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -5 -pf y
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +5 -pf y
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
	#grep "^CONECT" $e3 > ${e3%.*}_ligand2.pdb
	#cat ${e3%.*}_ligand2.pdb >> ${e3%.*}_ligand.pdb
        #pdb4amber -i ${e3%.*}_ligand.pdb -o ${e3%.*}_ligand_oprot_prep1.pdb --reduce
        #antechamber -fi pdb -i ${e3%.*}_ligand_oprot_prep1.pdb -fo mol2 -o ${e3%.*}_ligand_oprot_prep1_gas.mol2 -c gas -pf y
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc 0 -pf y -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -1 -pf y -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +1 -pf y -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -2 -pf y -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +2 -pf y -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -3 -pf y -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +3 -pf y -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +4 -pf y -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -4 -pf y -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -5 -pf y -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +5 -pf y -dr n
else
        :
fi


if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas.mol2 ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas.mol2 ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 0
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 1
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc -1
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 2
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc -2
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 3
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc -3
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas.mol2 ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas.mol2 ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 0 -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 1 -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc -1 -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 2 -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc -2 -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 3 -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc -3 -dr n
else
        :
fi

#sed -i '/X.*nan\|nan.*X/d' ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi
#sed -i '/X.*nan\|nan.*X/d' ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi

parmchk2 -f prepi -i ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
parmchk2 -f prepi -i ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep11.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@



tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m2 = loadpdb ${e3%.*}_r_prep11.pdb        # load receptor with ion and cofactor
check m2
charge m2
saveamberparm m2 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@

cp ${poi%.*}_r_prep11.pdb ${poi%.*}_r_prep00012.pdb
cp ${e3%.*}_r_prep11.pdb ${e3%.*}_r_prep00012.pdb


declare -i count=1
declare -i pcount=12
declare -i p2count=13

for i in $(seq 5);
do
if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Paramter expansion susbtring extraction ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis000${count}.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis000${count}.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis000${count}.out);
A="$(cut -d' ' -f1 <<<$del)"
a=${A:1:3}
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${poi%.*}_r_prep009.pdb >> ${poi%.*}_r_prep010.pdb
sed -e "/${a}.*${B}.*H/d" ${poi%.*}_r_prep000${pcount}.pdb >> ${poi%.*}_r_prep000${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep000${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :
fi;
done;

declare -i count=1
declare -i pcount=12
declare -i p2count=13

for i in $(seq 5);
do
if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Parameter expansion substring extraction ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis000${count}.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis000${count}.out);
C="$(cut -d' ' -f1 <<<$del2)"
c=${C:1:3}
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/H/d'  ${e3%.*}_r_prep009.pdb >> ${e3%.*}_r_prep010.pdb
sed -e "/${c}.*${D}.*H/d" ${e3%.*}_r_prep000${pcount}.pdb >> ${e3%.*}_r_prep000${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep000${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :
fi;
done;

cp ${poi%.*}_r_prep00017.pdb ${poi%.*}_r_prep1112.pdb
cp ${e3%.*}_r_prep00017.pdb ${e3%.*}_r_prep1112.pdb

declare -i pcount=12
declare -i p2count=13

for i in $(seq 35);
do

if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Hyd removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis11${count}.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis11${count}.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis11${count}.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

sed -e "/${A}.*${B}.*H/d" ${poi%.*}_r_prep11${pcount}.pdb >> ${poi%.*}_r_prep11${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep11${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :

fi;
done;

declare -i count=1
declare -i pcount=12
declare -i p2count=13

for i in $(seq 35);
do

if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Hyd removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis11${count}.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis11${count}.out);
C="$(cut -d' ' -f1 <<<$del2)"
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/H/d'  ${e3%.*}_r_prep006.pdb >> ${e3%.*}_r_prep007.pdb
sed -e "/${C}.*${D}.*H/d" ${e3%.*}_r_prep11${pcount}.pdb >> ${e3%.*}_r_prep11${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep11${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :

fi;
done;

cp ${poi%.*}_r_prep1145.pdb ${poi%.*}_r_prep11112.pdb
cp ${e3%.*}_r_prep1145.pdb ${e3%.*}_r_prep11112.pdb

declare -i count=1
declare -i pcount=12
declare -i p2count=13

for i in $(seq 5);
do
if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Paramter expansion susbtring extraction ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis111${count}.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis111${count}.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis111${count}.out);
A="$(cut -d' ' -f1 <<<$del)"
a=${A:1:3}
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${poi%.*}_r_prep009.pdb >> ${poi%.*}_r_prep010.pdb
sed -e "/${a}.*${B}.*H/d" ${poi%.*}_r_prep111${pcount}.pdb >> ${poi%.*}_r_prep111${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep111${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :
fi;
done;

declare -i count=1
declare -i pcount=12
declare -i p2count=13

for i in $(seq 5);
do
if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Parameter expansion substring extraction ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis111${count}.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis111${count}.out);
C="$(cut -d' ' -f1 <<<$del2)"
c=${C:1:3}
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/H/d'  ${e3%.*}_r_prep009.pdb >> ${e3%.*}_r_prep010.pdb
sed -e "/${c}.*${D}.*H/d" ${e3%.*}_r_prep111${pcount}.pdb >> ${e3%.*}_r_prep111${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep111${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :
fi;
done;

cp ${poi%.*}_r_prep11117.pdb ${poi%.*}_r_prep0012.pdb
cp ${e3%.*}_r_prep11117.pdb ${e3%.*}_r_prep0012.pdb

declare -i count=1
declare -i pcount=12
declare -i p2count=13

for i in $(seq 35);
do
if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - BACKBONE removal in progress ...";
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis00${count}.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis00${count}.out);
del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis00${count}.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"
C="$(cut -d' ' -f1 <<<$del2)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${poi%.*}_r_prep00${pcount}.pdb >> ${poi%.*}_r_prep00${p2count}.pdb
sed -e "/${C}.*${A}.*${B}/d" ${poi%.*}_r_prep00${pcount}.pdb >> ${poi%.*}_r_prep00${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep00${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :
fi;
done;

declare -i count=1
declare -i pcount=12
declare -i p2count=13

for i in $(seq 35);
do
if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - BACKBONE removal in progress ...";
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis00${count}.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis00${count}.out);
del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis00${count}.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"
C="$(cut -d' ' -f1 <<<$del2)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${e3%.*}_r_prep00${pcount}.pdb >> ${e3%.*}_r_prep00${p2count}.pdb
sed -e "/${C}.*${A}.*${B}/d" ${e3%.*}_r_prep00${pcount}.pdb >> ${e3%.*}_r_prep00${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m2 = loadpdb ${e3%.*}_r_prep00${p2count}.pdb        # load receptor with ion and cofactor
check m2
charge m2
saveamberparm m2 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :

fi;
done;



if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - O1P removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis001.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis001.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis001.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/P/d'  ${poi%.*}_r_prep0047.pdb >> ${poi%.*}_r_prep001.pdb
sed -e "/O1P.*${A}.*${B}/d" ${poi%.*}_r_prep0047.pdb >> ${poi%.*}_r_prep001.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep001.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - O1P removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis002.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis002.out);
C="$(cut -d' ' -f1 <<<$del2)"
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/P/d'  ${e3%.*}_r_prep0047.pdb >> ${e3%.*}_r_prep001.pdb
sed -e "/O1P.*${C}.*${D}/d" ${e3%.*}_r_prep0047.pdb >> ${e3%.*}_r_prep001.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep001.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - O2P removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis003.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis003.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis003.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/P/d'  ${poi%.*}_r_prep12.pdb >> ${poi%.*}_r_prep002.pdb
sed -e "/O2P.*${A}.*${B}/d" ${poi%.*}_r_prep001.pdb >> ${poi%.*}_r_prep002.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep002.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - O2P removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis004.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis004.out);
C="$(cut -d' ' -f1 <<<$del2)"
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/P/d'  ${e3%.*}_r_prep12.pdb >> ${e3%.*}_r_prep002.pdb
sed -e "/O2P.*${C}.*${D}/d" ${e3%.*}_r_prep001.pdb >> ${e3%.*}_r_prep002.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep002.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@
fi;


if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - O3P removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis005.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis005.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis005.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/P/d'  ${poi%.*}_r_prep12.pdb >> ${poi%.*}_r_prep003.pdb
sed -e "/O3P.*${A}.*${B}/d" ${poi%.*}_r_prep002.pdb >> ${poi%.*}_r_prep003.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep003.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - O3P removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis006.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis006.out);
C="$(cut -d' ' -f1 <<<$del2)"
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/P/d'  ${e3%.*}_r_prep12.pdb >> ${e3%.*}_r_prep003.pdb
sed -e "/O3P.*${C}.*${D}/d" ${e3%.*}_r_prep002.pdb >> ${e3%.*}_r_prep003.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep003.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@
fi;



if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - P removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis007.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis007.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis007.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/P/d'  ${poi%.*}_r_prep003.pdb >> ${poi%.*}_r_prep004.pdb
sed -e "/${A}.*${B}.*P/d" ${poi%.*}_r_prep003.pdb >> ${poi%.*}_r_prep004.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep004.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - P removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis008.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosi008.out);
C="$(cut -d' ' -f1 <<<$del2)"
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/P/d'  ${e3%.*}_r_prep003.pdb >> ${e3%.*}_r_prep004.pdb
sed -e "/${C}.*${D}.*P/d" ${e3%.*}_r_prep003.pdb >> ${e3%.*}_r_prep004.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep004.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - P removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis009.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis009.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis009.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/P/d'  ${poi%.*}_r_prep004.pdb >> ${poi%.*}_r_prep005.pdb
sed -e "/${A}.*${B}.*P/d" ${poi%.*}_r_prep004.pdb >> ${poi%.*}_r_prep005.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep005.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - P removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis010.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis010.out);
C="$(cut -d' ' -f1 <<<$del2)"
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/P/d'  ${e3%.*}_r_prep004.pdb >> ${e3%.*}_r_prep005.pdb
sed -e "/${C}.*${D}.*P/d" ${e3%.*}_r_prep004.pdb >> ${e3%.*}_r_prep005.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep005.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - P removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis011.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis011.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis011.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/P/d'  ${poi%.*}_r_prep005.pdb >> ${poi%.*}_r_prep006.pdb
sed -e "/${A}.*${B}.*P/d" ${poi%.*}_r_prep005.pdb >> ${poi%.*}_r_prep006.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep006.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - P removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis012.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis012.out);
C="$(cut -d' ' -f1 <<<$del2)"
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/P/d'  ${e3%.*}_r_prep005.pdb >> ${e3%.*}_r_prep006.pdb
sed -e "/${C}.*${D}.*P/d" ${e3%.*}_r_prep005.pdb >> ${e3%.*}_r_prep006.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep006.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Hyd removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis013.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis013.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis013.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${poi%.*}_r_prep006.pdb >> ${poi%.*}_r_prep007.pdb
sed -e "/${A}.*${B}.*H/d" ${poi%.*}_r_prep006.pdb >> ${poi%.*}_r_prep007.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep007.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Hyd removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis014.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis014.out);
C="$(cut -d' ' -f1 <<<$del2)"
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/H/d'  ${e3%.*}_r_prep006.pdb >> ${e3%.*}_r_prep007.pdb
sed -e "/${C}.*${D}.*H/d" ${e3%.*}_r_prep006.pdb >> ${e3%.*}_r_prep007.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep007.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Hyd removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis015.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis015.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis015.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${poi%.*}_r_prep007.pdb >> ${poi%.*}_r_prep008.pdb
sed -e "/${A}.*${B}.*H/d" ${poi%.*}_r_prep007.pdb >> ${poi%.*}_r_prep008.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep008.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Hyd removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis016.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis016.out);
C="$(cut -d' ' -f1 <<<$del2)"
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/H/d'  ${e3%.*}_r_prep007.pdb >> ${e3%.*}_r_prep008.pdb
sed -e "/${C}.*${D}.*H/d" ${e3%.*}_r_prep007.pdb >> ${e3%.*}_r_prep008.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep008.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Hyd removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis017.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis017.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis017.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${poi%.*}_r_prep008.pdb >> ${poi%.*}_r_prep009.pdb
sed -e "/${A}.*${B}.*H/d" ${poi%.*}_r_prep008.pdb >> ${poi%.*}_r_prep009.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep009.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Hyd removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis018.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis018.out);
C="$(cut -d' ' -f1 <<<$del2)"
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/H/d'  ${e3%.*}_r_prep008.pdb >> ${e3%.*}_r_prep009.pdb
sed -e "/${C}.*${D}.*H/d" ${e3%.*}_r_prep008.pdb >> ${e3%.*}_r_prep009.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep009.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Paramter expansion susbtring extraction ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis019.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis019.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis019.out);
A="$(cut -d' ' -f1 <<<$del)"
a=${A:1:3}
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${poi%.*}_r_prep009.pdb >> ${poi%.*}_r_prep010.pdb
sed -e "/${a}.*${B}.*H/d" ${poi%.*}_r_prep009.pdb >> ${poi%.*}_r_prep010.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep010.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Parameter expansion substring extraction ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis020.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis020.out);
C="$(cut -d' ' -f1 <<<$del2)"
c=${C:1:3}
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/H/d'  ${e3%.*}_r_prep009.pdb >> ${e3%.*}_r_prep010.pdb
sed -e "/${c}.*${D}.*H/d" ${e3%.*}_r_prep009.pdb >> ${e3%.*}_r_prep010.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep010.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - HIS removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis021.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis021.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis021.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${poi%.*}_r_prep010.pdb >> ${poi%.*}_r_prep011.pdb
sed -e "/HIS.*${B}.*H/d" ${poi%.*}_r_prep010.pdb >> ${poi%.*}_r_prep011.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep011.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - HIS removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis022.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis022.out);
C="$(cut -d' ' -f1 <<<$del2)"
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/H/d'  ${e3%.*}_r_prep010.pdb >> ${e3%.*}_r_prep011.pdb
sed -e "/HIS.*${D}.*H/d" ${e3%.*}_r_prep010.pdb >> ${e3%.*}_r_prep011.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep011.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@
fi;

declare -i count=23
declare -i pcount=11
declare -i p2count=12

for i in $(seq 35);
do
if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - BACKBONE removal in progress ...";
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis0${count}.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis0${count}.out);
del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis0${count}.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"
C="$(cut -d' ' -f1 <<<$del2)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${poi%.*}_r_prep0${pcount}.pdb >> ${poi%.*}_r_prep0${p2count}.pdb
sed -e "/${C}.*${A}.*${B}/d" ${poi%.*}_r_prep0${pcount}.pdb >> ${poi%.*}_r_prep0${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep0${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
	:
fi;
done;

declare -i count=24
declare -i pcount=11
declare -i p2count=12

for i in $(seq 35);
do
if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - BACKBONE removal in progress ...";
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis0${count}.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis0${count}.out);
del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis0${count}.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"
C="$(cut -d' ' -f1 <<<$del2)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${e3%.*}_r_prep0${pcount}.pdb >> ${e3%.*}_r_prep0${p2count}.pdb
sed -e "/${C}.*${A}.*${B}/d" ${e3%.*}_r_prep0${pcount}.pdb >> ${e3%.*}_r_prep0${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m2 = loadpdb ${e3%.*}_r_prep0${p2count}.pdb        # load receptor with ion and cofactor
check m2
charge m2
saveamberparm m2 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else 
	:

fi;
done;

if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - P, O1P, O2P, O3P removal in progress ..."

sed -i '/ATOM.*\<P\>/d' ${e3%.*}_r_prep046.pdb
sed -i '/ATOM.*\<O1P\>/d' ${e3%.*}_r_prep046.pdb
sed -i '/ATOM.*\<O2P\>/d' ${e3%.*}_r_prep046.pdb
sed -i '/ATOM.*\<O3P\>/d' ${e3%.*}_r_prep046.pdb



tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep046.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - P, O1P, O2P, O3P removal in progress ..."

sed -i '/ATOM.*\<P\>/d' ${poi%.*}_r_prep046.pdb
sed -i '/ATOM.*\<O1P\>/d' ${poi%.*}_r_prep046.pdb
sed -i '/ATOM.*\<O2P\>/d' ${poi%.*}_r_prep046.pdb
sed -i '/ATOM.*\<O3P\>/d' ${poi%.*}_r_prep046.pdb



tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep046.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@
fi;

declare -i count=1
declare -i pcount=46
declare -i p2count=47

for i in $(seq 8);
do
if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Paramter expansion susbtring extraction ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis22${count}.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis22${count}.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis22${count}.out);
A="$(cut -d' ' -f1 <<<$del)"
a=${A:1:3}
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${poi%.*}_r_prep009.pdb >> ${poi%.*}_r_prep010.pdb
sed -e "/${a}.*${B}.*H/d" ${poi%.*}_r_prep0${pcount}.pdb >> ${poi%.*}_r_prep0${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep0${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :
fi;
done;

declare -i count=1
declare -i pcount=46
declare -i p2count=47

for i in $(seq 8);
do
if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Parameter expansion substring extraction ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis22${count}.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis22${count}.out);
C="$(cut -d' ' -f1 <<<$del2)"
c=${C:1:3}
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/H/d'  ${e3%.*}_r_prep009.pdb >> ${e3%.*}_r_prep010.pdb
sed -e "/${c}.*${D}.*H/d" ${e3%.*}_r_prep0${pcount}.pdb >> ${e3%.*}_r_prep0${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep0${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :
fi;
done;


if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_special2.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@

fi;

if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_special2.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@

fi;

ambpdb -p ${poi%.*}_r_prep2.prmtop -c ${poi%.*}_r_prep2.inpcrd -mol2 -sybyl > ${poi%.*}_r_prep3.mol2
ambpdb -p ${e3%.*}_r_prep2.prmtop -c ${e3%.*}_r_prep2.inpcrd -mol2 -sybyl > ${e3%.*}_r_prep3.mol2

echo "STARTING PROPOSE PP DOCKING"

#bash push_renum.sh ${poi%.*}_r_prep1.pdb ${poi%.*}_r_prep1_renum.txt ${poi%.*}_r_prep3.mol2 > ${poi%.*}_r_prep3_renum.mol2
#if [ ! -f ${poi%.*}_r_prep3_renum.mol2 ] || [ ! -s ${poi%.*}_r_prep3_renum.mol2 ]; then
#        echo "contournement de renum!!!"
#        awk '!($2="")' ${poi%.*}_r_prep1_renum.txt >> ${poi%.*}_f_prep1_renum.txt
#        bash push_renum.sh ${poi%.*}_r_prep1.pdb ${poi%.*}_f_prep1_renum.txt ${poi%.*}_r_prep3.mol2 > ${poi%.*}_r_prep3_renum.mol2
#fi
#
#if [ ! -f ${poi%.*}_r_prep3_renum.mol2 ] || [ ! -s ${poi%.*}_r_prep3_renum.mol2 ]; then
#        echo "error is persisting in push_renum.sh for $poi"
#fi

#bash push_renum.sh ${e3%.*}_r_prep1.pdb ${e3%.*}_r_prep1_renum.txt ${e3%.*}_r_prep3.mol2 > ${e3%.*}_r_prep3_renum.mol2
#if [ ! -f ${e3%.*}_r_prep3_renum.mol2 ] || [ ! -s ${e3%.*}_r_prep3_renum.mol2 ]; then
#        awk '!($2="")' ${e3%.*}_r_prep1_renum.txt >> ${e3%.*}_f_prep1_renum.txt
#        bash push_renum.sh ${e3%.*}_r_prep1.pdb ${poi%.*}_f_prep1_renum.txt ${e3%.*}_r_prep3.mol2 > ${e3%.*}_r_prep3_renum.mol2
#fi
#
#if [ ! -f ${e3%.*}_r_prep3_renum.mol2 ] || [ ! -s ${e3%.*}_r_prep3_renum.mol2 ]; then
#        echo "error is persisting in push_renum.sh for $e3"
#fi


sed -i '/HETATM.*ZN\|ZN.*HETATM/d' ${poi%.*}_r_prep11.pdb
sed -i '/HETATM.*ZN\|ZN.*HETATM/d' ${e3%.*}_r_prep11.pdb
sed -i '/HETATM.*NI\|NI.*HETATM/d' ${poi%.*}_r_prep11.pdb
sed -i '/HETATM.*NI\|NI.*HETATM/d' ${e3%.*}_r_prep11.pdb
sed -i '/HETATM.*MG\|MG.*HETATM/d' ${poi%.*}_r_prep11.pdb
sed -i '/HETATM.*MG\|MG.*HETATM/d' ${e3%.*}_r_prep11.pdb
sed -i '/HETATM.*CA\|CA.*HETATM/d' ${poi%.*}_r_prep11.pdb
sed -i '/HETATM.*CA\|CA.*HETATM/d' ${e3%.*}_r_prep11.pdb

python ${PROTACable}/PROTACable_stage_2/get_pocket.py ${poi%.*}_r_prep11.pdb ${poi%.*}_r_pocket.txt
python ${PROTACable}/PROTACable_stage_2/get_pocket.py ${e3%.*}_r_prep11.pdb ${e3%.*}_r_pocket.txt

sed -i '$ s/.$//' ${poi%.*}_r_pocket.txt
sed -i '$ s/.$//' ${e3%.*}_r_pocket.txt

sed -i '1s/^/:/' ${poi%.*}_r_pocket.txt
sed -i '1s/^/:/' ${e3%.*}_r_pocket.txt

pocket=$(head -n 1 ${poi%.*}_r_pocket.txt)
pocket2=$(head -n 1 ${e3%.*}_r_pocket.txt)

echo "$pocket"
echo "$pocket2"

ambmask -p ${poi%.*}_r_prep2.prmtop -c ${poi%.*}_r_prep2.inpcrd -prnlev 0 -out pdb -find "$pocket"| awk '/ATOM/{print $2}' > ${poi%.*}_r_prep3.hit
ambmask -p ${e3%.*}_r_prep2.prmtop -c ${e3%.*}_r_prep2.inpcrd -prnlev 0 -out pdb -find "$pocket2"| awk '/ATOM/{print $2}' > ${e3%.*}_r_prep3.hit


#BIN=${PROPOSEHOME?"undefined"}/bin

${PROPOSEHOME}/ProPOSE TARGET=${poi%.*}_r_prep3.mol2 \
          LIGAND=${e3%.*}_r_prep3.mol2 \
          OUTTAR=${poi%.*}_r_prep_pp.mol2 \
          OUTLIG=${e3%.*}_l_prep_pp.mol2 \
          TARGET_PRM=${poi%.*}_r_prep2.prmtop \
          LIGAND_PRM=${e3%.*}_r_prep2.prmtop \
	  LIGAND_RMS=${e3%.*}_r_prep3.mol2 \
	  LICENSE=${LICENSE_PATH} \
          NOUT=100 \
	  NUM_THREAD=32
#	  

rm *out *special*pdb *prep1*pdb *prep0*pdb *hit *frcmod *pocket*txt *inpcrd *special*txt *sslink *frcmod *prepi *gas.mol2 *renum*txt
	
if [ -f ${poi%.*}_r_prep_pp.mol2 ] || [ -s ${poi%.*}_r_prep_pp.mol2 ]; then
	rm *special* *hit *prep2* *prep1* *ligand* sqm.* ANTECHAMBER* ATOMTYPE* *log *noHET* *sed* *sslink* *txt *renum*mol2
	echo "SUCCESS"
else
	echo "FAILED"
fi



rm *scores pp_scores.csv
awk '/SCORES:/,/Saved/' slurm*out  >> trial.scores
tail -n +2 trial.scores >> trial2.scores
head -n -1 trial2.scores > temp.txt ; mv temp.txt trial2.scores
column -t trial2.scores >> trial3.scores
sed -e 's/\s\+/,/g' trial3.scores > pp_scores.csv
#sed -i '2,102s/^.//' pp_scores.csv
rm *scores


rm *filtered* *fix*
#for p in *pp*mol2; do sed -e '/P.*UNL/!b' -e 's/C.3/P.3/g' $p >> ${p%.*}_p1fixed.mol2;done

${PROTACable}/PROTACable_stage_2/python filter_mol2.py


sh ${PROTACable}/PROTACable_stage_2/moe.sh $poi $e3

