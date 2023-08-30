import numpy as np
import torch.nn as nn


BASE_1HOT = {
    'A': np.array([1, 0, 0, 0]),
    'C': np.array([0, 1, 0, 0]),
    'G': np.array([0, 0, 1, 0]),
    'T': np.array([0, 0, 0, 1]),
    'W': np.array([.5, 0, 0, .5]),
    'S': np.array([0, .5, .5, 0]),
    'M': np.array([.5, .5, 0, 0]),
    'K': np.array([0, 0, .5, .5]),
    'R': np.array([.5, 0, .5, 0]),
    'Y': np.array([0, .5, 0, .5]),
    'B': np.array([0, 1. / 3, 1. / 3, 1. / 3]),
    'D': np.array([1. / 3, 0, 1. / 3, 1. / 3]),
    'H': np.array([1. / 3, 1. / 3, 0, 1. / 3]),
    'V': np.array([1. / 3, 1. / 3, 1. / 3, 0]),
    'N': np.array([.25, .25, .25, .25])
}


def one_hot_encode_seq(seq):
  encoded = np.zeros(shape = (4, len(seq)), dtype = np.float32)
  for i, base in enumerate(seq):
    encoded[:, i] = BASE_1HOT[base]
  return encoded


def get_DeepSEA_model():
  return nn.Sequential(
    #nn.Conv1d(
    #  in_channels = 4,
    #  out_channels = 320,
    #  kernel_size = 8,
    #  stride = 1
    #),
    nn.Conv1d(4, 16, 8, 1),
    nn.Threshold(
      threshold = 0,
      value = 1e-06
    ),
    nn.MaxPool1d(
      kernel_size = 4,
      stride = 4
    ),
    nn.Dropout(p = 0.2),
    #nn.Conv1d(320, 480, 8, 1),
    nn.Conv1d(16, 16, 8, 1),
    nn.Threshold(0, 1e-06),
    nn.MaxPool1d(4, 4),
    nn.Dropout(0.2),
    #nn.Conv1d(480, 960, 8, 1),
    nn.Conv1d(16, 16, 8, 1),
    nn.Threshold(0, 1e-06),
    nn.Dropout(0.5),
    nn.Flatten(),
    #nn.Linear(
    #  in_features = 50880,
    #  out_features = 925
    #),
    nn.Linear(848, 50),
    nn.Linear(50, 1),
    nn.Sigmoid()
  )

def snp_mutate(seq):
  seq = seq.upper()
  seqs = [seq]
  seq_list = list(seq)
  bases = {'A', 'C', 'T', 'G'}
  for i, base in enumerate(seq_list):
    other_bases = bases - {base}
    for other_base in other_bases:
      new_seq_list = seq_list.copy()
      new_seq_list[i] = other_base
      seqs.append("".join(new_seq_list))
  return seqs
