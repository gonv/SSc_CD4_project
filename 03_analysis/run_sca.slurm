#!/bin/bash
#SBATCH --job-name=SCA
#SBATCH --partition=bch-compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --time=00:20:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=gonzalo.villanuevamartin@childrens.harvard.edu
#SBATCH --output=/temp_work/ch258010/logs/00_CD4/out/05_SCA_%A.out
#SBATCH --error=/temp_work/ch258010/logs/00_CD4/err/05_SCA_%A.err

# --- env ---
source ~/.bashrc
conda activate condaR441 

# --- job ---
echo "[$(date +%H:%M:%S)] Running SCA job"
Rscript /lab-share/IM-Gutierrez-e2/Public/gvillanueva/00_CD4/scripts/run_sca.R
echo "[$(date +%H:%M:%S)] Job completed."
