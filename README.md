# Transcriptomic signatures of liver tumour recurrence

This repository contains the scripts used for the analyses presented in the manuscript:

"Integrative machine learning of multiomic data identifies the tumor adjacent liver microenvironment as a stable predictor of hepatocellular carcinoma recurrence"

## Overview

This repository contains the analysis scripts used to investigate the biological origin and stability of recurrence-associated molecular signals in hepatocellular carcinoma (HCC).

The study systematically compares tumor and adjacent (non-tumorous) liver tissue across independent cohorts to determine where reproducible and clinically relevant recurrence signals reside.

Key analyses include:
- Comparative analysis of tumor vs adjacent tissue for recurrence prediction
- Cross-cohort validation (HCC-KI ↔ HCC-TCGA) to assess generalizability
- Multiomic integration (mRNA + miRNA)
- Identification of coherent biological programs associated with recurrence, including:
  * Increased ribosome biogenesis and RNA processing
  * Decreased immune-related signaling
- Evaluation of signal stability and transferability across datasets

The results demonstrate that recurrence-associated molecular signals are more stable and biologically coherent in the tumor-adjacent liver microenvironment than in tumor tissue.

## Versions

R version 4.4.0

## Data availability

TCGA data can be obtained from: https://portal.gdc.cancer.gov/

Karolinska cohort data are available upon reasonable request.
