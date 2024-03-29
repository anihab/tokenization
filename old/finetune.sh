#!/bin/bash
#SBATCH -N 1
#SBATCH -p GPU
#SBATCH -t 2:00:00
#SBATCH --job-name=tokenization
#SBATCH --gpus=v100-32:8
#SBATCH -o /jet/home/ahabib/TokenizationFinetune%j.outerror

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

SCRIPT_PATH=/ocean/projects/bio230026p/ahabib/SCRIPTS/tokenization

# Arguments
# Arguments
BACTERIA_INPUT="/ocean/projects/bio230026p/ahabib/RAW_DATA/BACTERIA"
PHAGE_INPUT="/ocean/projects/bio230026p/ahabib/RAW_DATA/PHAGE"
OUTPUT_DIR="/ocean/projects/bio230026p/ahabib/RAW_DATA/formatted"
METHOD=kmer
KMER=6
VOCAB_INPUT=None

# Output Directories
BACTERIA_OUTPUT=$OUTPUT_DIR/BACTERIA
PHAGE_OUTPUT=$OUTPUT_DIR/PHAGE

# Create output directories
mkdir $BACTERIA_OUTPUT
mkdir $PHAGE_OUTPUT

echo "TIME: Start: = `date +"%Y-%m-%d %T"`"
cd $SCRIPT_PATH
pwd
if [[ $BACTERIA_INPUT == *.txt ]]; then
    python3 tokenization_list.py --b $BACTERIA_INPUT \
                             --p $PHAGE_INPUT \
                             --b_out $BACTERIA_OUTPUT \
                             --p_out $PHAGE_OUTPUT \
                             --method $METHOD \
                             --k $KMER \
                             --vocab $VOCAB_DIR
else
    python3 tokenization.py --b $BACTERIA_INPUT \
                        --p $PHAGE_INPUT \
                        --b_out $BACTERIA_OUTPUT \
                        --p_out $PHAGE_OUTPUT \
                        --method $METHOD \
                        --k $KMER \
                        --vocab $VOCAB_DIR
fi

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
