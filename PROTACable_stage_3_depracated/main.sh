#!/bin/bash

##SBATCH --partition=<YOUR_GPU_NODE_NAME>
#SBATCH --job-name=protacs

# example: $PROTACable/PROTACable_stage_3/main.sh <path_to_ternaries_output_from_stage2>
source ~/.bashrc
conda activate GNN_DOVE

GNN="${PROTACable}/PROTACable_stage_3/GNN_DOVE"

init_dir=$1

desired_dir="${init_dir%ternaries*}"

cd "$desired_dir" || exit 1 

for filepath in ternaries/*/top_20_pooled*; do
    filepath="$(pwd)/${filepath}"
    (cd $GNN && sh run_gnn.sh $filepath)
done;

input_dir="${PROTACable}/PROTACable_stage_3/GNN_DOVE/Predict_Result/Multi_Target/Fold_-1_Result/"
output_dir="ternaries/"
#process_dir="process_4/"


# Loop over all directories inside the input directory
for dir in ${input_dir}/top_20_pooled*/ ; do
    target=$(basename dir)
    subdir_name=$(basename "$dir")
    suffix="${subdir_name#top_20_pooled}"
    echo "$suffix"
    mkdir ${desired_dir}/${output_dir}/${suffix}/top_predicted_pose

    target=$(basename dir)
    # Check if the Predict_sort.txt file exists
    if [ -f "${dir}Predict_sort.txt" ]; then

        # Read the first line of the file and split it into an array
        read -a array <<< "$(head -n 1 ${dir}Predict_sort.txt)"

        #identifier1=${array[6]}_${array[7]}
        #full_identifier2=${array[10]}
        #identifier2=${full_identifier2%%-*}
        #echo "$identifier1 , $identifier2"
        identifier=$array

        # Copy files from the process directory to the output directory
        cp ${desired_dir}/${output_dir}/${suffix}/${target}/${identifier} ${desired_dir}/${output_dir}/${suffix}/top_predicted_pose

    fi
done


