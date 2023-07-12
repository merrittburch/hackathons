import logging
from pathlib import Path
import sys

from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import torch
import torch.nn as nn
import torch.optim as optim

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
training_intervals_path = Path(sys.argv[2])
param_file_path = Path(sys.argv[3])

if not reference_path.exists() or not training_intervals_path.exists():
  logging.error("Either the provided reference or training BED file doesn't exist")
  sys.exit(1)

# load our reference genome into a dictionary
ref = {rec.id: rec for rec in SeqIO.parse(reference_path, 'fasta')}

# read in the training data, one-hot encode the sequences
logging.info('One-hot encoding sequences')
samples = []

with open(training_intervals_path) as training_intervals_file:
  for line in training_intervals_file:
    chrom, start, end, name, bound = line.rstrip().split("\t")
    samples.append({
      '1hot_seq': one_hot_encode_seq(ref[chrom][int(start):int(end)]),
      'bound': int(bound)
    })

# take first 10% of training data as validation set
n_val = round(len(samples) / 10)
validation_set = samples[:n_val]
training_set = samples[n_val:]

#device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
device = 'cpu'

# load the model architecture and define our objective function and optimization algorithm
DeepSEA = get_DeepSEA_model().to(device)
criterion = nn.BCELoss()
optimizer = optim.SGD(DeepSEA.parameters(), lr = 0.1)

# begin our training loop
batch_size = 128
history = {
  'training': [],
  'validation': []
}
epochs = 3
for epoch in range(epochs):
  logging.info(f"Training epoch {epoch + 1}")
  running_training_loss = 0.0
  running_validation_loss = 0.0
  # split our training data into batches
  for batch_index_start in range(0, len(training_set), batch_size):
    batch = training_set[batch_index_start:batch_index_start + batch_size]
    DeepSEA = DeepSEA.train()

    training_inputs = torch.tensor([sample['1hot_seq'] for sample in batch], dtype=torch.float32, device = device)
    training_targets = torch.tensor([sample['bound'] for sample in batch], dtype=torch.float32, device = device).reshape(-1, 1)
    
    optimizer.zero_grad()

    training_outputs = DeepSEA.forward(training_inputs)
    training_loss = criterion(training_outputs, training_targets)
    training_loss.backward()
    optimizer.step()

    running_training_loss += training_loss.item()

  # compute validation loss every epoch
  for batch_index_start in range(0, len(validation_set), batch_size):
    batch = validation_set[batch_index_start:batch_index_start + batch_size]
    DeepSEA = DeepSEA.eval()
    with torch.no_grad():
      validation_inputs = torch.tensor([sample['1hot_seq'] for sample in batch], dtype=torch.float32)
      validation_targets = torch.tensor([sample['bound'] for sample in batch], dtype=torch.float32).reshape(-1, 1)
      validation_outputs = DeepSEA.forward(validation_inputs)
    running_validation_loss += criterion(validation_outputs, validation_targets).item()

  avg_training_loss = running_training_loss / len(training_set)
  avg_validation_loss = running_validation_loss / len(validation_set)
  history['training'].append(avg_training_loss)
  history['validation'].append(avg_validation_loss)
  logging.info(f"Average training loss: {avg_training_loss:.4f}")
  logging.info(f"Average validation loss: {avg_validation_loss:.4f}")

# save model in TorchScript format
DeepSEA.eval()
DeepSEA_traced = torch.jit.trace(DeepSEA, torch.randn(1, 4, 1000))
DeepSEA_traced.save('ABF1_DeepSEA.pt')

# plot training and validation loss curves
fig, ax = plt.subplots()
ax.plot(range(1, epochs + 1), history['training'], label = 'Training', color = 'orange')
ax.plot(range(1, epochs + 1), history['validation'], label= 'Validation', color = 'blue')
ax.set_xlabel('Epoch')
ax.set_ylabel('Loss')
ax.legend()
ax.set_title('Training Curves')
ax.get_xaxis().set_major_locator(MaxNLocator(integer = True))
fig.savefig('training.png')

# save model to file
torch.save(DeepSEA.state_dict(), param_file_path)

