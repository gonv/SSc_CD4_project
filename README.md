# Transcriptional landscape of CD4+ T cells in Systemic Sclerosis

**Authors:** Gonzalo Villanueva-Martin, Gonzalo Borrego-Yaniz, Marialbert Acosta-Herrera,  José Luis Callejas-Rubio, Norberto Ortego, Norbert Mages, Stefan Börno, Maria Gutierrez-Arcelus, Javier Martin,* and Lara Bossini-Castillo*

*These authors contributed equally.

**Institution:** Instituto de Parasitologia y Biomedicina, Lopez-Neyra / Boston Children's Hospital / Harvard Medical School

**Repository author**: Gonzalo Villanueva-Martin

---

## 🧬 Project overview

This repository contains the full single-cell RNA sequencing (scRNA-seq) analysis pipeline for the study:

> **“Transcriptional landscape of CD4+ T cells in Systemic Sclerosis”**

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
```
---
## 🖥️ Shiny App
This resource is available as a Shiny App [docs-link]: https://mgalab.shinyapps.io/CD4inSScResource/

---
## 📚 References
<a id="1">[1]</a> 
Villanueva-Martin G, et al. (2026). Transcriptional landscape of CD4+ T cells in Systemic Sclerosis. bioRxiv.

