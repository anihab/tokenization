#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=soc-gpu-np
#SBATCH --partition=soc-gpu-np
#SBATCH --job-name=dnabert2
#SBATCH --time=1:00:00
#SBATCH -o /uufs/chpc.utah.edu/common/home/u1049062/OUT

set -x

# Load Modules
nvidia-smi 
module load cuda
nvidia-smi
module load miniconda3

which python3
python --version
which pip
conda list

# echo "starting dna env on conda"
# conda activate dna
# which python3
# python --version
# which pip
# conda list

data_path="/uufs/chpc.utah.edu/common/home/u1323098/anisa/RAW_DATA/SPLIT"
output_path="/uufs/chpc.utah.edu/common/home/u1323098/anisa/TOKENIZED_DATA/DNABERT2"

echo "The provided data_path is $data_path"
echo "TIME: Start: = `date +"%Y-%m-%d %T"`"

python3 /uufs/chpc.utah.edu/common/home/u1323098/anisa/MODELS/DNABERT_2/finetune/train.py \
            --model_name_or_path zhihan1996/DNABERT-2-117M \
            --data_path  $data_path \
            --kmer -1 \
            --run_name DNABERT2_TEST \
            --model_max_length 128 \
            --per_device_train_batch_size 8 \
            --per_device_eval_batch_size 16 \
            --gradient_accumulation_steps 1 \
            --learning_rate 2e-4 \
            --num_train_epochs 3 \
            --fp16 \
            --save_steps 200 \
            --output_dir $output_path \
            --evaluation_strategy steps \
            --eval_steps 200 \
            --warmup_steps 50 \
            --logging_steps 100000 \
            --overwrite_output_dir True \
            --log_level info \
            --find_unused_parameters False

echo "TIME: End: = `date +"%Y-%m-%d %T"`" 