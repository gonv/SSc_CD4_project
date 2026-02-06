#!/bin/bash
#SBATCH --job-name=pseudobulk_clusters
#SBATCH --partition=bch-compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=6-10
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=00:30:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=gonzalo.villanuevamartin@childrens.harvard.edu
#SBATCH --output=/temp_work/ch258010/logs/00_CD4/out/10_pseudobulk_ct_%j_%A_%a.out
#SBATCH --error=/temp_work/ch258010/logs/00_CD4/err/10_pseudobulk_ct_%j_%A_%a.err

# --- env ---
source ~/.bashrc
conda activate condaR441

# --- script options ---
#/lab-share/IM-Gutierrez-e2/Public/gvillanueva/00_CD4/scripts/run_pseudobulk_clusters.R - # for the clusters log 09 #SBATCH --array=0-11
#/lab-share/IM-Gutierrez-e2/Public/gvillanueva/00_CD4/scripts/run_pseudobulk_cellTypes.R - # for the cell types log 10 #SBATCH --array=1-5

# --- vars ---
#CON="$SLURM_ARRAY_TASK_ID"
CT_FILE="/lab-share/IM-Gutierrez-e2/Public/gvillanueva/00_CD4/scripts/cellTypes.txt"
CELLTYPE=$(awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR == ARRID { print $1 }' "$CT_FILE")
echo "Array task ID: $SLURM_ARRAY_TASK_ID"
echo "Cell type seleccionado: $CELLTYPE"

# --- run --
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running pseudobulk job for cell type: $CELLTYPE"
Rscript /lab-share/IM-Gutierrez-e2/Public/gvillanueva/00_CD4/scripts/run_pseudobulk_cellTypes.R "$CELLTYPE"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Job completed."