#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=soc-gpu-np
#SBATCH --partition=soc-gpu-np
#SBATCH --job-name=tokenization
#SBATCH --time=1:00:00
#SBATCH -o /uufs/chpc.utah.edu/common/home/u1323098/TokenizationTest%j.outerror

# Load Modules 
module purge
module use $HOME/MyModules
module load miniconda3/latest
echo "starting TOKENIZE env on conda"
source activate TOKENIZE
which python3
python --version 
which pip
conda list

SCRIPT_PATH="/uufs/chpc.utah.edu/common/home/u1323098/anisa/SCRIPTS"

# Arguments
BACTERIA_DIR="/uufs/chpc.utah.edu/common/home/u1323098/anisa/RAW_DATA/BACTERIA/FASTA"
PHAGE_DIR="/uufs/chpc.utah.edu/common/home/u1323098/anisa/RAW_DATA/PHAGE/FASTA"
OUTPUT_DIR="/uufs/chpc.utah.edu/common/home/u1323098/anisa/TOKENIZED_DATA/BPE"
METHOD="bpe"
K=None
VOCAB=None

echo "TIME: Start: = `date +"%Y-%m-%d %T"`"
cd $SCRIPT_PATH
pwd
python3 tokenization.py --b $BACTERIA_DIR --p $PHAGE_DIR --o1 $OUTPUT_DIR --method $METHOD --vocab $VOCAB
echo "TIME: End: = `date +"%Y-%m-%d %T"`"

# notes:
# test sample size - 1.2MB - 5 bacteria files, 5 phage files
# codon -- 1 to 5 sec
# 6mer -- 1 to 5 sec
# bpe -- 8 min
#     -- overhead pre-train full dataset = 10 sec 