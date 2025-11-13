#!/bin/bash
#SBATCH --output=log/%j.out
#SBATCH --job-name=ele
#SBATCH --get-user-env=L
#SBATCH --exclusive
##SBATCH --partition bulk
#SBATCH --time=24:00:00
#SBATCH --ntasks=32
source ~/.bashrc
conda activate capstone-melting
python ele.py "$@"
