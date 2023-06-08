#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail

# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-06-05
# Updated... 2023-06-05
#
# Description:
# Create conda environment, download hugging face transformer model, set up 
# jupyter notebook within environment
# ------------------------------------------------------------------------------

# Create a conda environment
conda create -n hugging-face -c conda-forge transformers

# Enter environment with a name
conda activate hugging-face

# Install the model from the hugging face source library
pip install --upgrade git+https://github.com/huggingface/transformers.git

# Check to see if pytorch is using the gpu
python3
import torch
torch.cuda.is_available() # it works if returns true

# Check to see if we're running from the correct conda environment
which jupyter-lab
jupyter-kernelspec list

# Install the ipykernel
python -m ipykernel install --user --name hugging-face --display-name "hugging-face"

# Start a screen environment
screen -S jupyter

# Start a jupyter notebook
jupyter-lab --ip=0.0.0.0 --port=8017 --no-browser
jupyter-lab --ip=0.0.0.0 --port=8017

# Copy the link
# Exit out of the screen session 