# Libraries
import os
import gzip
import argparse
import random

import pandas as pd
from Bio import SeqIO

"""\
Given directories of bacteria and phage fasta files, parse and format each file into a
.csv that includes the sequence of given max_length and label (0 for bacteria, 1 for phage)

Arguments:
  bacteria_dir -- str, path to directory of bacteria fasta files
  phage_dir -- str, path to directory of phage fasta files
  max_length -- int, max sequence length for parsing
  shift_amount -- int, length n to shift sequences by
  random_start -- bool, whether to parse with random start positions instead of 0
  b_out -- str, path to directory for all or bacteria output
  p_out -- str, path to directory for phage ouput
"""
def read_files(bacteria_dir, phage_dir, max_length, shift_amount, random_start, **kwargs):
  # read and verify arguments
  b_out = kwargs.get('b_out', None)
  p_out = kwargs.get('p_out', None)

  if p_out is None:
    p_out = b_out

  # format all files in bacteria_dir 
  for filename in os.listdir(bacteria_dir):
    f = os.path.join(bacteria_dir, filename)
    if os.path.isfile(f) and not is_file_read(b_out, filename): # TODO: fix to skip not fasta
      format(f, max_length, shift_amount, random_start, b_out)     
  # format all files in phage_dir
    for filename in os.listdir(phage_dir):
      f = os.path.join(phage_dir, filename)
      if os.path.isfile(f) and not is_file_read(p_out, filename): # TODO: fix to skip not fasta
        format(f, max_length, shift_amount, random_start, p_out)
  
"""\
Runs fasta files through tokenizer and adds the label of 1 for phage and
0 for bacteria. Then shuffles the rows in the dataframe and saves to CSV 

Arguments:
  filepath -- str, path to fasta file
  max_length -- int, max sequence length for parsing
  output_dir -- str, path to directory for output
"""
def format(filepath, max_length, shift_amount, random_start, output_dir):
  filename = os.path.basename(filepath)
  filename = filename.split('.')[0]

  # process data to get sequences of appropriate length and save output to csv
  if shift_amount is None:
    df = preprocess_data(filepath, max_length, random_start)
    write_csv(filename, df, output_dir)

"""\
Read fasta file and truncate sequences to appropriate length, returns dataframe

Arguments:
  filepath -- str, path to fasta file
  max_length -- int, maximum sequence length
Returns:
  df -- dataframe, includes the file title, > input line, start position, and sequence
""" 
def preprocess_data(filepath, max_length, random_start): 
  records = []

  f = filepath
  if filepath.endswith('.gz'):
    f = gzip.open(filepath, 'rt', encoding='utf-8')

  try:
    for record in SeqIO.parse(f, 'fasta'):
      filename = os.path.basename(filepath)
      name = filename.split('.')[0]
      segment = str(record.name)
      seq = str(record.seq).upper()
      pos = 0

      # if random_start, start from a random position in the sequence
      if random_start is True:
        r = random.randint(0, len(seq) - max_length)
        seq = seq[r:]

      # truncate sequences if longer than max_length
      if max_length is not None:
        while len(seq) > max_length:
          records.append(                  # add subsequence up to max_length
            {
              'name': name,
              'segment': segment,
              'start': pos,
              'sequence': seq[:max_length]
            }
          )
          seq = seq[max_length:]           # sequence continuing from max_length
          pos += max_length
      records.append(
          {
            'name': name,
            'segment': segment,
            'start': pos,
            'sequence': seq
          }
      )
  finally:
    df = pd.DataFrame(data=records)
    return df
  
"""\
Save the given dataframe to a _formatted csv file, which includes every
sequence and corresponding label.

Arguments:
  filename -- str, name of file being tokenized
  df -- dataframe, full dataframe of tokenized sequences
  output_dir -- str, path to directory for output
"""
def write_csv(filename, df, output_dir):
  formatted = df[['sequence', 'label']]
  formatted.to_csv(output_dir + '/' + filename + '_formatted.csv', encoding='utf-8', index=False, header=False) 

'''\
Determines whether or not a file has already been processed by checking
if the output filename exists in the output directory and has a size
greater than 0.

Arguments:
  filename -- str, the full name of the file to check
Returns:
  isfile -- bool, whether the file has already been tokenized
'''
def is_file_read(directory, filename):
  file_path = os.path.join(directory, filename.split('.')[0]  + '_tokenized.csv')
  if os.path.isfile(file_path) and os.path.getsize(file_path) > 0:
    return True
  else:
    return False

## MAIN -----------------------------------------------------------------------------------------------------  

def main():
  parser = argparse.ArgumentParser()
  # Parameters
  parser.add_argument(
        "--b", default=None, type=str, required=True, help="The input bacteria directory."
    )
  parser.add_argument(
        "--p", default=None, type=str, required=True, help="The input phage directory."
    )
  parser.add_argument(
        "--max_length", default=None, type=int, required=True, help="The max sequence length for parsing"
  )
  parser.add_argument(
        "--shift_amount", default=None, type=int, required=True, help="The amount n to shift sequences by when parsing, if 0 then will not shift"
  )
  parser.add_argument(
        "--random_start", default=None, type=bool, required=True, help="If true, the start locations of sequences will be randomly selected"
  )
  parser.add_argument(
        "--b_out", default=None, type=str, required=False, help="The first output directory, for bacteria if using both."
    )
  parser.add_argument(
        "--p_out", default=None, type=str, required=False, help="The second output directory, for phage if using both."
    )
  args = parser.parse_args()

  # read and format files
  read_files(bacteria_dir=args.b, 
             phage_dir=args.p, 
             max_length=args.max_length,
             shift_amount=args.shift_amount,
             random_start=args.random_start,
             b_out=args.b_out, 
             p_out=args.p_out)

if __name__ == "__main__":
    main()