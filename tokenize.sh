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

SCRIPT_PATH=/ocean/projects/bio230026p/ahabib/SCRIPTS/tokenization

# Arguments
BACTERIA_INPUT=
PHAGE_INPUT=
OUTPUT_DIR=
METHOD=
KMER=
VOCAB_INPUT=

echo "TIME: Start: = `date +"%Y-%m-%d %T"`"
cd $SCRIPT_PATH
pwd
if [[ $BACTERIA_INPUT == *.txt ]]; then
    python3 tokenization_list.py --b $BACTERIA_INPUT \
                                 --p $PHAGE_INPUT \
                                 --b_out $OUTPUT_DIR \
                                 --method $METHOD \
                                 --vocab $VOCAB_INPUT
else
    python3 tokenization.py --b $BACTERIA_INPUT \
                            --p $PHAGE_INPUT \
                            --b_out $OUTPUT_DIR \
                            --method $METHOD \
                            --vocab $VOCAB_INPUT
fi
echo "TIME: End: = `date +"%Y-%m-%d %T"`" 
