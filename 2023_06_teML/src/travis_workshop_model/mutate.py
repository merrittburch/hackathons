import logging
from pathlib import Path
import sys

from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import torch

from lib import one_hot_encode_seq, get_DeepSEA_model, snp_mutate

logging.basicConfig(
  level=logging.INFO,
  format="{asctime} [{levelname}] {message}",
  style="{"
)
if len(sys.argv) != 3:
  logging.error("Wrong number of arguments provided")
  sys.exit(1)

seq_path = Path(sys.argv[1])
param_file_path = Path(sys.argv[2])

if not seq_path.exists() or not param_file_path.exists():
  logging.error("Either the provided sequence or param file doesn't exist")
  sys.exit(1)

# load in our sequence to be mutated
seq = str(SeqIO.read(seq_path, 'fasta').seq)

logging.info("Mutating sequence")
# mutate the sequence with all possible single SNPs
seqs = snp_mutate(seq)

logging.info("One-hot encoding sequences")
# one-hot encode the sequences
seqs = torch.tensor([one_hot_encode_seq(seq) for seq in seqs], dtype = torch.float32)

# load our pre-trained model
DeepSEA = get_DeepSEA_model()
DeepSEA.load_state_dict(torch.load(param_file_path))

# run our original and mutated sequences through the model in batches
DeepSEA = DeepSEA.eval()
with torch.no_grad():
  unmutated = DeepSEA.forward(seqs[0].reshape(1, 4, -1)).flatten()

logging.info("Predicting mutated sequences")
batch_size = 128
mutated = []
for batch_index_start in range(1, len(seqs), batch_size):
  batch = seqs[batch_index_start:batch_index_start + batch_size]
  DeepSEA = DeepSEA.eval()
  with torch.no_grad():
    inputs = batch
    outputs = DeepSEA.forward(inputs)
  mutated.append(outputs)

mutated = torch.cat(mutated).flatten()

logging.info("Calculating variant effects")
# compute variant effects
variant_effects = mutated - unmutated

logging.info("Plotting variant effects")
# plot the variant effects
fig, ax = plt.subplots()

for i in range(0, len(variant_effects), 3):
  position = i // 3
  plt.scatter([position] * 3, variant_effects[i:i+3])

# plot motif locations
for motif_loc in [846, 313, 641]:
  plt.axvline(motif_loc, linestyle="--", color = "red")

plt.title("in silico mutagenesis")
plt.xlabel("Position")
plt.ylabel("Predicted Binding Effect")
fig.tight_layout()
fig.savefig('mutagenesis.png')
