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

BACTERIA_DIR="/uufs/chpc.utah.edu/common/home/u1323098/anisa/RAW_DATA/BACTERIA/FASTA"
PHAGE_DIR="/uufs/chpc.utah.edu/common/home/u1323098/anisa/RAW_DATA/PHAGE/FASTA"
OUTPUT_DIR="/uufs/chpc.utah.edu/common/home/u1323098/anisa/TOKENIZED_DATA/BPE"

METHOD="bpe"
K=None

echo "TIME: Start: = `date +"%Y-%m-%d %T"`"
cd $SCRIPT_PATH
pwd
python3 tokenization.py --b $BACTERIA_DIR --p $PHAGE_DIR --o1 $OUTPUT_DIR --method $METHOD --k $K
echo "TIME: End: = `date +"%Y-%m-%d %T"`"

# notes:
# codon -- 6 sec for 5 bactera, 5 phage
# 6mer -- 6 sec for 5 bactera, 5 phage
# bpe -- 8.01 min for 5 bactera, 5 phage