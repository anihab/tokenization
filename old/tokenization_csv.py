# Libraries
import argparse
import gzip
import os

import pandas as pd
import numpy as np

from Bio import SeqIO
from tokenizers import Tokenizer, models, trainers, normalizers
from transformers import AutoModel, AutoTokenizer

# Globals
MAX_TOKENS = 510
VOCAB_SIZE = 50000

## READ INPUT -----------------------------------------------------------------------------------------------------  

"""\
Given a bacteria directory and a phage directory, tokenize all files
according to method of choice.

Arguments:
  input_path -- str, path to directory/list or a single csv file of selected, formatted DNA segments
  method -- str, tokenization method of choice
  k -- int, length of k if using kmer tokenization
  vocab -- str, path to directory of fasta files to build model vocab
  output_path -- str, path to directory for output csv
"""
def read_files(input_path, output_path, method, **kwargs):
  # read and verify arguments
  k = kwargs.get('k', None)
  vocab = kwargs.get('vocab', None)

  if method == 'kmer' and k == None:
    print("Missing argument, kmer tokenization method requires parameter \'k\'.")
  if method == 'bpe' and vocab == None:
    print("Missing argument, bpe tokenization method requires parameter \'vocab\'.")

  # build model vocabulary if using bpe
  if method == 'bpe':
    build_vocab(vocab)

  # tokenize all of the sequences
  if input_path.endswith('.txt'):
    read_list(input_path, method, k, output_path)
  elif input_path.endswith('.csv'):
    read_csv(input_path, method, k, output_path)
  else:
    read_directory(input_path, method, k, output_path)

"""\
Read all of the files in a given .txt list of files, and tokenize accordingly.
"""
def read_list(input_list, method, k, output_path):
  if os.path.isfile(input_list): 
    with open(input_list, 'r') as list:
      for f in list:
        f = f.strip()
        filename = os.path.basename(f) 
        if os.path.isfile(f) and not is_file_read(output_path, filename):
          tokenize(f, method, k, output_path)

"""\
Read all of the files in a given directory of files, and tokenize accordingly.
"""
def read_directory(input_dir, method, k, output_path):
  for filename in os.listdir(input_dir):
    f = os.path.join(input_dir, filename)
    if os.path.isfile(f) and not is_file_read(output_path, filename):
      tokenize(f, method, k, output_path)

"""\
Read all of the segments in a given csv file, and tokenize accordingly.
"""
def read_csv(input_csv, method, k, output_path):
  filename = os.path.basename(input_csv)
  if os.path.isfile(input_csv) and not is_file_read(output_path, filename):
      tokenize(input_csv, method, k, output_path)

'''\
Determines whether or not a file has already been processed by checking
if the output filename exists in the output directory and has a size
greater than 0.

Arguments:
  filename -- str, the full name of the file to check
Returns:
  isfile -- bool, whether the file has already been tokenized
'''
def is_file_read(output_path, filename):
  filepath = os.path.join(output_path, filename.split('.')[0]  + '_tokenized.csv')
  if os.path.isfile(filepath) and os.path.getsize(filepath) > 0:
    return True
  else:
    return False

"""\
Runs fasta files through tokenizer and adds the label of 1 for phage and
0 for bacteria. Then shuffles the rows in the dataframe and saves to CSV 

Arguments:
  filepath -- str, path to fasta file
  label -- int, 0 for bacteria or 1 for phage
  method -- str, tokenization method of choice
  output_dir -- str, path to directory for output
  k -- int, length of k if using kmer tokenization
"""
def tokenize(filepath, method, k, output_dir):
  sequences = []
  tokens = []

  filename = os.path.basename(filepath)
  filename = filename.split('.')[0]

  df = pd.read_csv(filepath)
  sequences = df['sequence'].values.tolist()

  # tokenize according to chosen method
  for seq in range(len(sequences)):
    if method == 'codon':
      tokens.append(seq2codon(sequences[seq]))
    elif method == 'kmer':
      tokens.append(seq2kmer(sequences[seq], k))
    elif method == 'bpe':
      tokens.append(seq2bpe(sequences[seq]))
  df['tokenized'] = tokens
  
  # shuffle and save to csv
  df = df.sample(frac=1).reset_index(drop=True)
  write_csv(filename, df, output_dir)
  return df

"""\
Save the given dataframe to two separate csv files:
1. _full.csv includes the name, start position, sequence, tokenized
   sequence, and label.
2. _tokenized.csv includes the tokenized sequence and the label.

Arguments:
  filename -- str, name of file being tokenized
  df -- dataframe, full dataframe of tokenized sequences
  output_dir -- str, path to directory for output
""" 
def write_csv(filename, df, output_dir):
  # df.to_csv(directory + "/" + filename + '_full.csv', encoding='utf-8', index=False)
  tokenized = df[['tokenized', 'label']]
  tokenized.to_csv(output_dir + '/' + filename + '_tokenized.csv', encoding='utf-8', index=False, header=False, sep='\t')

## TOKENIZATION METHODS ----------------------------------------------------------------------------------------------------- 

"""\
Convert a sequence to codons

Arguments:
  seq -- str, original sequence
Returns:
  codons -- str, codons separated by space
"""
def seq2codon(seq):
  codon = [seq[i:i+3] for i in range(0,len(seq),3)]
  codons = " ".join(codon)
  return codons

"""\
Convert a sequence to kmers

Arguments:
  seq -- str, original sequence
  k -- int, kmer of length k
Returns:
  kmers -- str, kmers separated by space
"""
def seq2kmer(seq, k):
  kmer = [seq[i:i+k] for i in range(len(seq)+1-k)]
  kmers = " ".join(kmer)
  return kmers

"""\
Convert a sequence to byte pair encodings

Arguments:
  seq -- str, original sequence
Returns:
  output -- str, decoded tokens separated by a space
"""
def seq2bpe(sequence):                         
  #tokenizer = PreTrainedTokenizerFast(tokenizer_file="tokenizer.json")r
  tokenizer = AutoTokenizer.from_pretrained("tokenizer.json")
  encoded_input = tokenizer(sequence, return_tensors="pt")
  token_ids = encoded_input.input_ids
  output = " ".join(tokenizer.batch_decode(token_ids)) 
  return output

"""\
Build the vocabulary for the BPE model

Arguments:
  vocab -- str, directory or list of files to build model vocabulary on
"""
def build_vocab(vocab):
  sequences = parse_vocab(vocab)
  tokenizer = Tokenizer(models.BPE())
  # Customize the tokenizer to handle DNA sequences
  tokenizer.normalizer = normalizers.Sequence([normalizers.NFKC()])
  # Train the tokenizer on your DNA sequences
  trainer = trainers.BpeTrainer(vocab_size=VOCAB_SIZE)
  tokenizer.train_from_iterator(sequences, trainer=trainer)
  tokenizer.save("tokenizer.json")

"""\
Parse vocabulary for the BPE model to work with

Arguments:
  vocab -- str, directory or list of files to build model vocabulary on
Returns:
  sequences -- list, sequences to train model on
"""
def parse_vocab(vocab):
  sequences = []
  # If the input vocabulary is a directory
  if os.path.isdir(vocab):
    for filename in os.listdir(vocab):
      f = os.path.join(vocab, filename)
      if os.path.isfile(f):
        if f.endswith('.gz'):
          f = gzip.open(f, 'rt', encoding='utf-8')
        for record in SeqIO.parse(f, 'fasta'):
          sequences.append(str(record.seq).upper())
  # If the input vocabulary is a list of files
  elif os.path.isfile(vocab):
    with open(vocab, 'r') as list:
      for f in list:
        f = f.strip()
        filename = os.path.basename(f)
        if os.path.isfile(f):
          if f.endswith('.gz'):
            f = gzip.open(f, 'rt', encoding='utf-8')
          for record in SeqIO.parse(f, 'fasta'):
            sequences.append(str(record.seq).upper())
  return sequences

## MAIN -----------------------------------------------------------------------------------------------------  

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument(
        "--input", default=None, type=str, required=True, help="A path to directory/list or a single csv file of selected, formatted bacteria DNA segments"
    )
  parser.add_argument(
        "--output", default=None, type=str, required=False, help="The output directory path."
    )
  parser.add_argument(
        "--method", default=None, type=str, required=True, help="The tokenization method of choice: kmer, codon, or bpe."
    )
  parser.add_argument(
        "--k", default=None, type=int, required=False, help="Length k for kmer tokenization."
    )
  parser.add_argument(
        "--vocab", default=None, type=str, required=False, help="The directory or list of files to build the model vocabulary for if using the bpe method."
    )
  args = parser.parse_args()

  # read and tokenize files
  read_files(input_path=args.input,
             output_path=args.output,
             method=args.method, 
             k=args.k, 
             vocab=args.vocab)

if __name__ == "__main__":
    main()
