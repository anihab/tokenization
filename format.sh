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

# Locations
# NOTE: This script assumes you have an output directory set-up with BACTERIA and PHAGE sub-folders
#       OUTPUT_DIR ls 
#           BACTERIA
#           PHAGE
BACTERIA="/uufs/chpc.utah.edu/common/home/u1323098/anisa/RAW_DATA/BACTERIA"
PHAGE="/uufs/chpc.utah.edu/common/home/u1323098/anisa/RAW_DATA/PHAGE"
OUTPUT="/uufs/chpc.utah.edu/common/home/sundar-group2/ANISA/TOKENIZED_DATA/test"

FORMAT_LENGTH=500       # maximum sequence length for when parsing
SELECT_AMOUNT=12000     # number of bacteria and phage samples to respectively select for final csv

echo "TIME: Start: = `date +"%Y-%m-%d %T"`"
echo "begin formatting individual files"

cd $SCRIPT_PATH
if [[ $BACTERIA == *.txt ]]; then
    python3 tokenization_list.py --b $BACTERIA \
                        --p $PHAGE \
                        --b_out $OUTPUT/BACTERIA \
                        --p_out $OUTPUT/PHAGE \
                        --method format \
                        --format_length $FORMAT_LENGTH
else
    python3 tokenization.py --b $BACTERIA \
                        --p $PHAGE \
                        --b_out $OUTPUT/BACTERIA \
                        --p_out $OUTPUT/PHAGE \
                        --method format \
                        --format_length $FORMAT_LENGTH
fi                      

echo "completed formatting individual files, begin combining and shuffling data"

# Concat and shuffle the data
cd $OUTPUT/BACTERIA
cat *.csv > bacteria_unshuf.csv
shuf bacteria_unshuf.csv > bacteria_shuf.csv

cd $OUTPUT/PHAGE
cat *.csv > phage_unshuf.csv
shuf phage_unshuf.csv > phage_shuf.csv 

# Select a number n (n=SELECT_AMOUNT) of the sequences
head -n $SELECT_AMOUNT $OUTPUT/BACTERIA/bacteria_shuf.csv > $OUTPUT/BACTERIA/bacteria_selected.csv
head -n $SELECT_AMOUNT $OUTPUT/PHAGE/phage_shuf.csv > $OUTPUT/PHAGE/phage_selected.csv

echo "TIME: End: = `date +"%Y-%m-%d %T"`"