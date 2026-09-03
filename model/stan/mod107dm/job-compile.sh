#!/bin/bash -l

#SBATCH --job-name=cmdstanr_conda
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2gb
#SBATCH --time=00:05:00
#SBATCH --array=1-1%1
#SBATCH -p test

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# load
module load trytonp/r/4.5.2

export MODEL_SRC_DIR="/users/project1/pt01268/cmdstanr"

cd ${MODEL_SRC_DIR}

Rscript compile_stan.R