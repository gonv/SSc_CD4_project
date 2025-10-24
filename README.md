# SSc_CD4_project
Single-cell RNA-seq analysis of CD4⁺ T cells in Systemic Sclerosis

**Author:** Gonzalo Villanueva-Martin  
**Institution:** Boston Children's Hospital / Harvard Medical School  
**Project status:** Ongoing (Private repository until publication)

---

## 🧬 Project overview

This repository contains the full single-cell RNA sequencing (scRNA-seq) analysis pipeline for the study:

> **“Single-Cell RNA Sequencing Reveals Impaired Regulatory T Cell Function and a Pro-inflammatory Th17 Profile in Systemic Sclerosis”**

The project focuses on characterizing the transcriptional programs of **CD4⁺ T cells** in **Systemic Sclerosis (SSc)** compared to healthy controls, using **Seurat (R)** and **Scanpy (Python)** frameworks.  
All steps — from preprocessing to integration, clustering, and pathway analysis — are included.

---

## ⚙️ Environment

- **R 4.4.1**
- **Seurat 5.0+**
- **sctransform**, **harmony**, **fgsea**, **tidyverse**
- **Python 3.10+**
- **scanpy**, **scikit-learn**, **pandas**, **scipy**
- Cluster: **BCH Compute (SLURM)**

To reproduce the R environment:

```bash
conda create -n condaR441 r-base=4.4.1 r-seurat r-tidyverse r-harmony r-sctransform r-fgsea
