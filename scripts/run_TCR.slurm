#!/bin/bash
#SBATCH --job-name=CD4_TCR
#SBATCH --partition=bch-compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=80G
#SBATCH --time=12:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=gonzalo.villanuevamartin@childrens.harvard.edu
#SBATCH --output=/temp_work/ch258010/logs/00_CD4/out/01_TCR_%j.out
#SBATCH --error=/temp_work/ch258010/logs/00_CD4/err/01_TCR_%j.err

# --- env ---
source ~/.bashrc
conda activate condaScanpy

# Prevent BLAS oversubscription
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

python /lab-share/IM-Gutierrez-e2/Public/gvillanueva/00_CD4/scripts/run_TCR.py