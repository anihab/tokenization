#!/bin/bash
#SBATCH -N 1
#SBATCH -p GPU
#SBATCH -t 2:00:00
#SBATCH --job-name=split
#SBATCH --gpus=v100-32:8

set -x

module load anaconda3
echo "starting dna env on conda"
conda activate dna
which python3
python --version
which pip
conda list

# Set the input file path
input_file="/ocean/projects/bio230026p/ahabib/RAW_DATA/formatted.csv"

echo "TIME: Start: = `date +"%Y-%m-%d %T"`"

# Set the output file paths for train, test, and dev splits
train_file="/ocean/projects/bio230026p/ahabib/RAW_DATA/train.csv"
test_file="/ocean/projects/bio230026p/ahabib/RAW_DATA/test.csv"
dev_file="/ocean/projects/bio230026p/ahabib/RAW_DATA/dev.csv"

# Set the train-test-dev split ratios
train_ratio=0.7
test_ratio=0.15
dev_ratio=0.15

# Run Python script to perform train-test-dev split
python_script=$(cat << END
import pandas as pd
from sklearn.model_selection import train_test_split

# Load the input CSV file
data = pd.read_csv("$input_file")

# Perform train-test-dev split
train_data, temp_data = train_test_split(data, train_size=$train_ratio, test_size=$test_ratio+ $dev_ratio)

# Perform test-dev split on the remaining data
test_data, dev_data = train_test_split(temp_data, train_size=$test_ratio/($test_ratio + $dev_ratio), test_size=$dev_ratio/($test_ratio + $dev_ratio))

# Save train, test, and dev data to separate CSV files
train_data.to_csv("$train_file", index=False)
test_data.to_csv("$test_file", index=False)
dev_data.to_csv("$dev_file", index=False)
END
)

# Execute the Python script
python -c "$python_script"

echo "TIME: End: = `date +"%Y-%m-%d %T"`" 