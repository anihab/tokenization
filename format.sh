#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=soc-gpu-np
#SBATCH --partition=soc-gpu-np
#SBATCH --job-name=tokenization
#SBATCH --time=1:00:00
#SBATCH -o /uufs/chpc.utah.edu/common/home/u1323098/Format%j.outerror

# Load Modules 
module purge
module load anaconda3
echo "starting dna env on conda"
conda activate dna
which python3
python --version 
which pip
conda list

BACTERIA=/uufs/chpc.utah.edu/common/home/sundar-group2/ANISA/RAW_DATA/BACTERIA/formatted
PHAGE=/uufs/chpc.utah.edu/common/home/sundar-group2/ANISA/RAW_DATA/PHAGE/formatted
OUTPUT=/uufs/chpc.utah.edu/common/home/sundar-group2/ANISA/RAW_DATA

echo "TIME: Start: = `date +"%Y-%m-%d %T"`"

# Concat and shuffle the data
cd $BACTERIA
cat *.csv > bacteria_unshuf.csv
shuf bacteria_unshuf.csv > bacteria_shuf.csv

cd $PHAGE
cat *.csv > phage_unshuf.csv
shuf phage_unshuf.csv > phage_shuf.csv 

# Select a number n (n=30000 for now) of the sequences
head -n 30000 $BACTERIA/bacteria_shuf.csv
head -n 30000 $PHAGE/phage_shuf.csv

# Concat the phage_train.tsv and the bacteria_train.tsv together and shuffle again
cd $OUTPUT
cat  $BACTERIA/bacteria_shuf.csv $PHAGE/phage_shuf.csv > formatted_unshuf.csv
shuf formatted_unshuf.csv > formatted.csv

echo "TIME: End: = `date +"%Y-%m-%d %T"`"