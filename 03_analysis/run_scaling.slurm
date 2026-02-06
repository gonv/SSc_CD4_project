#!/bin/bash
#SBATCH --job-name=CD4_adata
#SBATCH --partition=bch-compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --time=12:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=gonzalo.villanuevamartin@childrens.harvard.edu
#SBATCH --output=/temp_work/ch258010/logs/00_CD4/out/03_scaling_%A.out
#SBATCH --error=/temp_work/ch258010/logs/00_CD4/err/03_scaling_%A.err

# --- env ---
source ~/.bashrc
conda activate condaR441 

echo "[$(date +%H:%M:%S)] Running scaling job"
Rscript /lab-share/IM-Gutierrez-e2/Public/gvillanueva/00_CD4/scripts/run_scaleing.R
echo "[$(date +%H:%M:%S)] Job completed."
