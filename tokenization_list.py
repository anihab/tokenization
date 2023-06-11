# -*- coding: utf-8 -*-

# Libraries
import argparse
import gzip
import os

import pandas as pd
import numpy as np

from Bio import SeqIO
from tokenizers import Tokenizer, models, trainers, normalizers
from transformers import PreTrainedTokenizerFast

# Globals
MAX_TOKENS = 510
BACTERIA_OUTPUT = ""
PHAGE_OUTPUT = ""

## Tokenize sequences given a list of files as input

"""\
Given lists of bacteria and phage fasta files, tokenize all files
according to method of choice.

Input:
  bacteria_files -- str, list of bacteria fasta files
  phage_files -- str, list of phage fasta files
  method -- str, tokenization method of choice
  k -- int, length of k if using kmer tokenization
  vocab -- str, path to directory of fasta files to build model vocab
"""
def read_files(bacteria_files, phage_files, method, **kwargs):
  k = kwargs.get('k', None)
  vocab = kwargs.get('vocab', None)
  # Ensure that input is correct
  if method == 'kmer' and k == None:
    print("Missing argument, kmer tokenization method requires parameter \'k\'.")
  if method == 'bpe' and vocab == None:
    print("Missing argument, bpe tokenization method requires parameter \'vocab\'.")
  # Build model vocabulary if using byte tokenization
  if method == 'bpe':
    build_vocab(vocab)
  # Tokenize all files in bacteria_files
  for f in bacteria_files:
    if os.path.isfile(f):
      tokenize(f, 0, method, k=k)
  # Tokenize all files in phage_files
  for f in phage_files:
    if os.path.isfile(f):
        tokenize(f, 1, method, k=k)

"""\
Runs fasta files through tokenizer and adds the label of 1 for phage and
0 for bacteria. Then shuffles the rows in the dataframe and saves to CSV 

Input:
  filepath -- str, path to fasta file
  label -- int, 0 for bacteria or 1 for phage
  method -- str, tokenization method of choice
  k -- int, length of k if using kmer tokenization
"""
def tokenize(filepath, label, method, **kwargs):
  sequences = []
  tokens = []

  k = kwargs.get('k', None)
  filename = os.path.basename(filepath)
  filename = filename.split('.')[0]

  # Calculate max segment length
  if method == 'codon':
    max_length = MAX_TOKENS * 3
  elif method == 'kmer':
    max_length = MAX_TOKENS - (k - 1)
  elif method == 'bpe':
    max_length = MAX_TOKENS

  # Process data to get sequences of appropriate length 
  df = preprocess_data(filepath, max_length)
  sequences = df['sequence'].values.tolist()

  # Tokenize according to chosen method
  for seq in range(len(sequences)):
    if method == 'codon':
      tokens.append(seq2codon(sequences[seq]))
    elif method == 'kmer':
      tokens.append(seq2kmer(sequences[seq], k))
    elif method == 'bpe':
      tokens.append(seq2bpe(sequences[seq]))
  df['tokenized'] = tokens
  df['label'] = [label] * len(tokens)
  
  # Shuffle and save to csv
  df = df.sample(frac=1).reset_index(drop=True)
  write_csv(filename, label, df)
  return df

"""\
Read fasta file and truncate sequences to appropriate length, returns dataframe

Input:
  filepath -- str, path to fasta file
  max_length -- int, maximum sequence length

Returns:
  df -- dataframe, includes the file title, > input line, start position, and sequence
""" 
def preprocess_data(filepath, max_length): 
  records = []

  f = filepath
  if filepath.endswith('.gz'):
    f = gzip.open(filepath, 'rt', encoding='utf-8')

  for record in SeqIO.parse(f, 'fasta'):
    filename = os.path.basename(filepath)
    name = filename.split('.')[0]
    segment = str(record.name)
    seq = str(record.seq).upper()
    pos = 0 
    # Truncate sequences if longer than max_length
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
  df = pd.DataFrame(data=records)
  return df

"""\
Read in sequences and tokens to attach labels and return dataframe

Input:
  sequences -- list, original sequences
  tokens -- list, tokenized sequences
  label -- int, 1 for phage or 0 for bacteria

Returns:
  df -- dataframe
""" 
def attach_labels(sequences, tokens, label):
  d = []
  for i in range(len(tokens)):
    d.append(
        {
          'sequence': sequences[i],
          'tokenized': tokens[i],
          'label': label
        }
    )
  df = pd.DataFrame(data=d)
  return df

"""\
Save the given dataframe to two separate csv files:
1. _full.csv includes the name, start position, sequence, tokenized
   sequence, and label.
2. _tokenized.csv includes the tokenized sequence and the label.

Input:
  df -- dataframe, full dataframe of tokenized sequences
"""
def write_csv(filename, label, df):
  if label == 0:
    directory = BACTERIA_OUTPUT
  else:
    directory = PHAGE_OUTPUT
  # df.to_csv(directory + "/" + filename + '_full.csv', encoding='utf-8', index=False)
  tokenized = df[['tokenized', 'label']]
  tokenized.to_csv(directory + "/" + filename + '_tokenized.csv', encoding='utf-8', index=False, header=False, sep='\t')
    

## Different tokenization methods

"""\
Convert a sequence to codons

Input:
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

Input:
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

Input:
  seq -- str, original sequence

Returns:
  output -- str, decoded tokens separated by a space
"""
def seq2bpe(sequence):                         
  tokenizer = PreTrainedTokenizerFast(tokenizer_file="dna_tokenizer.json")
  encoded_input = tokenizer(sequence, return_tensors="pt")
  token_ids = encoded_input.input_ids
  output = " ".join(tokenizer.batch_decode(token_ids)) 
  return output

"""\
Build vocabulary for the byte tokenization model to work with

Input:
  vocab -- str, directory or list of files to build model vocabulary on
"""
def build_vocab(vocab):
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
  else:
    for filename in vocab:
      if os.path.isfile(f):
        if f.endswith('.gz'):
          f = gzip.open(f, 'rt', encoding='utf-8')
        for record in SeqIO.parse(f, 'fasta'):
          sequences.append(str(record.seq).upper())
  train_bpe_tokenizer(sequences)

"""\
'Pretrain' or build the vocabulary for the BPE model

Input:
  sequences -- list, sequences to train model on
"""
def train_bpe_tokenizer(sequences):
  tokenizer = Tokenizer(models.BPE())

  # Customize the tokenizer to handle DNA sequences
  tokenizer.normalizer = normalizers.Sequence([normalizers.NFKC()])

  # Train the tokenizer on your DNA sequences
  trainer = trainers.BpeTrainer(vocab_size=50000)
  tokenizer.train_from_iterator(sequences, trainer=trainer)

  tokenizer.save("dna_tokenizer.json")

## MAIN

def main():
  parser = argparse.ArgumentParser()
  # Parameters
  parser.add_argument(
        "--b", default=None, type=str, required=True, help="The input bacteria files."
    )
  parser.add_argument(
        "--p", default=None, type=str, required=True, help="The input phage files."
    )
  parser.add_argument(
        "--o1", default=None, type=str, required=True, help="The first output directory, for bacteria if using both."
    )
  parser.add_argument(
        "--o2", default=None, type=str, required=False, help="The second output directory, for phage if using both."
    )
  parser.add_argument(
        "--method", default=None, type=str, required=True, help="The tokenization method of choice: kmer, codon, or bpe."
    )
  parser.add_argument(
        "--k", default=None, type=int, required=False, help="Length k for kmer tokenization."
    )
  parser.add_argument(
        "--vocab", default=None, type=str, required=False, help="The directory or list of files to build the model vocabulary, if using bpe."
    )
  args = parser.parse_args()

  global BACTERIA_OUTPUT
  BACTERIA_OUTPUT = args.o1

  global PHAGE_OUTPUT
  if args.o2 != None: 
    PHAGE_OUTPUT = args.o2
  else:
    PHAGE_OUTPUT = args.o1

  read_files(bacteria_files=args.b, phage_files=args.p, method=args.method, k=args.k, vocab=args.vocab)

if __name__ == "__main__":
    main()