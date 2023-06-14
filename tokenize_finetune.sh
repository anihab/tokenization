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
BACTERIA_DIR="/uufs/chpc.utah.edu/common/home/u1323098/sundar-group-space2/PHAGE/DATASETS/BACTERIA_RAW/FASTA/ncbi-genomes-2023-03-09/TEST"
PHAGE_DIR="/uufs/chpc.utah.edu/common/home/u1323098/anisa/RAW_DATA/PHAGE/FASTA"
METHOD="kmer"
KMER=6
VOCAB_DIR=None

# Output Directories
BACTERIA_OUTPUT=$BACTERIA_DIR/TOKENIZED
PHAGE_OUTPUT=$PHAGE_DIR/TOKENIZED

# Create output directories
mkdir "$BACTERIA_OUTPUT"
mkdir "$PHAGE_OUTPUT"

echo "TIME: Start: = `date +"%Y-%m-%d %T"`"
cd $SCRIPT_PATH
pwd
python3 tokenization.py --b $BACTERIA_DIR --p $PHAGE_DIR --o1 $BACTERIA_OUTPUT --o2 $PHAGE_OUTPUT --method $METHOD --k $KMER --vocab $VOCAB_DIR

# Concat all of the files together into one file
cd $OUTPUT_DIR
pwd
cd $BACTERIA_OUTPUT
cat *.csv > bacteria_train_unshuffled.csv
cd $PHAGE_OUTPUT
cat *.csv > phage_train_unshuffled.csv 

# Shuffle the data
cd $BACTERIA_OUTPUT
shuf -o bacteria_train.tsv bacteria_train_unshuffled.csv
cd $PHAGE_OUTPUT
shuf -o phage_train.tsv phage_train_unshuffled.csv

# Select a number n (n=15000 for now) of the sequences
head -n 15000 $BACTERIA_OUTPUT/bacteria_train.tsv
head -n 15000 $PHAGE_OUTPUT/phage_train.tsv

# Concat the phage_train.tsv and the bacteria_train.tsv together and shuffle again
cd $OUTPUT_DIR
cat  $BACTERIA_OUTPUT/bacteria_train.tsv $PHAGE_OUTPUT/phage_train.tsv > train_unshuffled.tsv
shuf train_unshuffled.tsv > train_shuffled.tsv

# Add the header (sequence        label) tab in between
echo -e "sequence\tlabel" | cat - train_shuffled.tsv > train.tsv

echo "TIME: End: = `date +"%Y-%m-%d %T"`"