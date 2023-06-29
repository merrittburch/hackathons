#!/bin/bash
#SBATCH --job-name=23-candida
#SBATCH --account=PAS2510
#SBATCH --time=10:00:00
#SBATCH --nodes=1 --ntasks-per-node=28
#SBATCH --mail-user=mbb262@cornell.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL