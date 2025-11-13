#!/bin/bash
#SBATCH --output=time.out
#SBATCH --job-name=pi_8
#SBATCH --get-user-env=L
#SBATCH --exclusive
##SBATCH --partition bulk
#SBATCH --time=24:00:00
#SBATCH --ntasks=128
source ~/.bashrc
conda activate capstone-melting
calphy_kernel -i input.yaml -k $1


