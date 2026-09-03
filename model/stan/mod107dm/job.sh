#!/bin/bash -l

#SBATCH --job-name=cmdstan_conda_batch
#SBATCH -N 1
#SBATCH --cpus-per-task=48
#SBATCH -p batch
#SBATCH --time=00:05:00
#SBATCH --mem=24gb
#SBATCH --array=1-10%10

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# load
module load trytonp/r/4.5.2

export MODEL_SRC_DIR="/users/project1/pt01268/cmdstanr"

cd ${MODEL_SRC_DIR}

Rscript run_evsi.R

