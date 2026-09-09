# Genomic Characterization of Stealth Methicillin-Resistant *Staphylococcus aureus* (OS-MRSA) — Analysis Code & Minimal Dataset

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue.svg)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-%E2%89%A53.9-blue.svg)](https://www.python.org/)

This repository contains the reproducible analysis code, figure generation scripts, and de-identified minimal dataset supporting the study:

> **"Genomic characterization of stealth methicillin-resistant *Staphylococcus aureus* (OS-MRSA) from companion animal wound infections in Dhaka, Bangladesh — A One Health concern"**  
> Target Journal: *One Health Advances* (BioMed Central / Springer Nature)

---

## 🔬 Study Summary

Standard phenotypic susceptibility tests frequently fail to detect oxacillin-susceptible *mecA*-positive *Staphylococcus aureus* (OS-MRSA; "stealth" MRSA). In this cross-sectional hospital survey of 107 companion animals (cats, dogs, horse) presenting with clinical wound infections in Dhaka, Bangladesh:
* **Prevalence:** *S. aureus* was confirmed in 22.4% (24/107) and *mecA*-positive MRSA in 16.8% (18/107) of animals.
* **Stealth Phenotype:** 33.3% (6/18) of *mecA*-positive isolates were phenotypically susceptible to cefoxitin by disc diffusion.
* **Genomic Architecture:** Short-read whole-genome sequencing (WGS) of the index feline OS-MRSA isolate confirmed a class B *mec* complex (*mecA*-$\Delta$*mecR1*, lacking *mecI*) co-assembled with *ccrA*/*ccrB* homologs on a 73-kb contig (NODE_13), with human-associated immune evasion cluster genes (IEC type D: *scn*, *sak*, *sea*).

---

## 📁 Repository Structure

```
OS-MRSA-Analysis/
├── data/
│   ├── Additional_file_1_Deidentified_Minimal_Dataset.csv   # Primary de-identified animal & isolate dataset (n = 107)
│   └── Additional_file_1_Deidentified_Minimal_Dataset.xlsx  # Excel workbook version with variable definitions
├── scripts/
│   ├── reproduce_table4_fisher_tests.R                      # Base R script reproducing Table 4 inferential statistics
│   └── create_figures.py                                   # Python script reproducing publication figures (Figs 1–4, S1–S4)
├── README.md                                               # Project overview and reproduction instructions
└── LICENSE                                                 # MIT License
```

---

## 📊 Reproducing Table 4 (Inferential Statistics in R)

The exact conditional odds ratios (OR), 95% confidence intervals, and two-sided Fisher's exact test $p$-values reported in **Table 4** can be reproduced using **Base R** (no external packages required):

```bash
# Clone the repository
git clone https://github.com/nehalhasnain/OS-MRSA-Analysis.git
cd OS-MRSA-Analysis

# Run the R analysis script
Rscript scripts/reproduce_table4_fisher_tests.R
```

This reproduces all 9 subgroup comparisons bit-for-bit and exports `Table_4_reproduced.csv`.

---

## 📈 Reproducing Publication Figures (Python)

Figures 1–4 and Supplementary Figures S1–S4 are generated using Python:

```bash
# Install dependencies
pip install matplotlib pandas numpy scipy

# Run figure generation script
python scripts/create_figures.py
```

Outputs are saved at 300/600 dpi TIFF, PNG, and vector PDF formats.

---

## 🧬 Public Genomic Data

All raw sequencing reads and assembled genomic data are deposited in public repositories:
* **NCBI BioProject:** [PRJNA1328371](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1328371)
* **NCBI BioSample:** [SAMN51332128](https://www.ncbi.nlm.nih.gov/biosample/SAMN51332128)
* **NCBI Sequence Read Archive (SRA):** [SRR35391355](https://www.ncbi.nlm.nih.gov/sra/SRR35391355)
* **GenBank Whole Genome Shotgun Assembly:** [DBKGHL000000000.1](https://www.ncbi.nlm.nih.gov/nuccore/DBKGHL000000000.1)

---

## ⚖️ License

This project is licensed under the [MIT License](LICENSE).
