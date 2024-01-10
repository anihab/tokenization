#!/bin/bash
#SBATCH -N 1
#SBATCH -p GPU
#SBATCH -t 2:00:00
#SBATCH --job-name=tokenization
#SBATCH --gpus=v100-32:8
#SBATCH -o /jet/home/ahabib/Tokenization%j.outerror

# Load Modules 
module purge
module load anaconda3
echo "starting dna env on conda"
conda activate dna
which python3
python --version 
which pip
conda list

SCRIPT_PATH="/ocean/projects/bio230026p/ahabib/SCRIPTS/tokenization"

# Arguments
BACTERIA_INPUT="/ocean/projects/bio230026p/ahabib/RAW_DATA/BACTERIA"
PHAGE_INPUT="/ocean/projects/bio230026p/ahabib/RAW_DATA/PHAGE"
OUTPUT_DIR1="/ocean/projects/bio230026p/ahabib/RAW_DATA/BACTERIA/formatted"
OUTPUT_DIR2="/ocean/projects/bio230026p/ahabib/RAW_DATA/PHAGE/formatted"
METHOD=format
KMER=0
VOCAB_INPUT=None

echo "TIME: Start: = `date +"%Y-%m-%d %T"`"

cd $SCRIPT_PATH
pwd
if [[ $BACTERIA_INPUT == *.txt ]]; then
    python3 tokenization_list.py --b $BACTERIA_INPUT \
                                 --p $PHAGE_INPUT \
                                 --b_out $OUTPUT_DIR1 \
                                 --p_out $OUTPUT_DIR2 \
                                 --method $METHOD \
                                 --k $KMER \
                                 --vocab $VOCAB_INPUT
else
    python3 tokenization.py --b $BACTERIA_INPUT \
                            --p $PHAGE_INPUT \
                            --b_out $OUTPUT_DIR1 \
                            --p_out $OUTPUT_DIR2 \
                            --method $METHOD \
                            --k $KMER \
                            --vocab $VOCAB_INPUT
fi
echo "TIME: End: = `date +"%Y-%m-%d %T"`" 
