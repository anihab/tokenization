#!/bin/bash

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

SCRIPT_PATH="/uufs/chpc.utah.edu/common/home/u1323098/anisa/SCRIPTS/tokenization"

BACTERIA="/uufs/chpc.utah.edu/common/home/u1323098/anisa/RAW_DATA/SHIFT/BACTERIA"
PHAGE="/uufs/chpc.utah.edu/common/home/u1323098/anisa/RAW_DATA/SHIFT/PHAGE"
OUTPUT="/uufs/chpc.utah.edu/common/home/sundar-group2/ANISA/RAW_DATA/SHIFT"

MAX_LENGTH=10          # maximum sequence length for when parsing
SHIFT_AMOUNT=5
RANDOM=False
SELECT_AMOUNT=12000     # number of bacteria and phage samples to respectively select for final csv

echo "TIME: Start: = `date +"%Y-%m-%d %T"`"

cd $SCRIPT_PATH
echo "begin formatting individual files"

python3 pretokenization.py --b $BACTERIA \
                        --p $PHAGE \
                        --max_length $MAX_LENGTH \
                        --random_start $RANDOM \
                        --shift_amount $SHIFT_AMOUNT \
                        --b_out $OUTPUT/BACTERIA \
                        --p_out $OUTPUT/PHAGE
                                              
echo "completed formatting individual files, begin combining and shuffling data"

# Concat and shuffle the data for BACTERIA
cd $OUTPUT/BACTERIA
cat *.csv > bacteria_unshuf.csv
shuf bacteria_unshuf.csv > bacteria_shuf.csv
rm bacteria_unshuf.csv

# Concat and shuffle the data for PHAGE
cd $OUTPUT/PHAGE
cat *.csv > phage_unshuf.csv
shuf phage_unshuf.csv > phage_shuf.csv
rm phage_unshuf.csv

# Select a number n (n=SELECT_AMOUNT) of the sequences
head -n $SELECT_AMOUNT $OUTPUT/BACTERIA/bacteria_shuf.csv > $OUTPUT/BACTERIA/bacteria_selected.csv
head -n $SELECT_AMOUNT $OUTPUT/PHAGE/phage_shuf.csv > $OUTPUT/PHAGE/phage_selected.csv

echo "completed combining and shuffling data"

echo "TIME: End: = `date +"%Y-%m-%d %T"`"