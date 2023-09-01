import logging
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from sklearn.metrics import (
  average_precision_score,
  precision_recall_curve,
  roc_curve
)


logging.basicConfig(
  level=logging.INFO,
  format="{asctime} [{levelname}] {message}",
  style="{"
)

# check our command-line arguments
if len(sys.argv) != 3:
  logging.error("Wrong number of arguments provided")
  sys.exit(1)

test_intervals_path = Path(sys.argv[1])
predictions_path = Path(sys.argv[2])

if not test_intervals_path.exists() or not predictions_path.exists():
  logging.error("Either the test BED file or the predictions file doesn't exist")
  sys.exit(1)

# read in the test data
truth = []
with open(test_intervals_path) as test_intervals_file:
  for line in test_intervals_file:
    truth.append(int(line.rstrip().split("\t")[4]))

# read in the predictions
predictions = []
with open(predictions_path) as predictions_file:
  for line in predictions_file:
    predictions.append(float(line.rstrip()))

# Use scikit-learn to evaluate our classification results
logging.info(f"Average precision score: {average_precision_score(truth, predictions)}")

# plot a PR curve
fig, ax = plt.subplots()

## determine the baseline
n_true = truth.count(1)
base_y = n_true / len(truth)

precision, recall, thresholds = precision_recall_curve(truth, predictions)
ax.plot(recall, precision)
ax.plot([0, 1], [base_y, base_y], color = "gray", linestyle="--", label = "Baseline")
ax.set_xlim((0, 1))
ax.set_ylim((0, 1))
ax.set_xlabel("Recall")
ax.set_ylabel("Precision")
ax.legend()

fig.savefig('pr.png')

# plot a ROC curve
fig, ax = plt.subplots()

fpr, tpr, thesholds = roc_curve(truth, predictions)
ax.plot(fpr, tpr)
ax.plot([0, 1], [0, 1], color = "gray", linestyle="--", label="Baseline")
ax.set_xlim((0, 1))
ax.set_ylim((0, 1))
ax.set_xlabel("False Positive Rate")
ax.set_ylabel("True Positive Rate")
ax.legend()

fig.savefig('roc.png')

