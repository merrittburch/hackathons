import logging
from pathlib import Path
import sys

from Bio import SeqIO
import torch

from lib import one_hot_encode_seq, get_DeepSEA_model


logging.basicConfig(
  level=logging.INFO,
  format="{asctime} [{levelname}] {message}",
  style="{"
)

# check our command-line arguments
if len(sys.argv) != 4:
  logging.error("Wrong number of arguments provided")
  sys.exit(1)

reference_path = Path(sys.argv[1])
test_intervals_path = Path(sys.argv[2])
param_file_path = Path(sys.argv[3])

if not reference_path.exists() or not test_intervals_path.exists() or not param_file_path.exists():
  logging.error("One or more input files don't exist")
  sys.exit(1)

# load our reference genome
ref = {rec.id: rec for rec in SeqIO.parse(reference_path, 'fasta')}

# read in the test data, one-hot encode the sequences
logging.info("One-hot encoding sequences")
samples = []
with open(test_intervals_path) as test_intervals_file:
  for line in test_intervals_file:
    chrom, start, end, name, bound = line.rstrip().split("\t")
    samples.append({
      '1hot_seq': one_hot_encode_seq(ref[chrom][int(start):int(end)]),
      'bound': int(bound)
    })

# load our pre-trained model
DeepSEA = get_DeepSEA_model()
DeepSEA.load_state_dict(torch.load(param_file_path))

# run our test data through the pre-trained model
# the procedure is similar to the validation set
batch_size = 128
for batch_index_start in range(0, len(samples), batch_size):
  batch = samples[batch_index_start:batch_index_start + batch_size]
  DeepSEA = DeepSEA.eval()
  with torch.no_grad():
    test_inputs = torch.tensor([sample['1hot_seq'] for sample in batch], dtype = torch.float32)
    test_outputs = DeepSEA.forward(test_inputs)
    sys.stdout.write("\n".join([str(float(e)) for e in test_outputs.flatten()]) + "\n")


