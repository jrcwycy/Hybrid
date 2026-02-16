# Hybrid

Repository with pipelines and utilities for gathering, integrating, and analyzing single-cell RNA sequencing data from transcription factor (TF) perturbation experiments. Datasets include one control condition (normal proliferating human fibroblasts) and three perturbation conditions: MYOD1 overexpression (OE), PRRX1 knockdown (KD), and sequential PRRX1 KD / MYOD1 OE. Current workflows cover preprocessing and quality control. The repository will be updated as the project is ongoing.


## Repository Structure

- `pipelines/` - Collection of workflow scripts
    - `EPI2ME/` - Nextflow pipeline for TF perturbation experiments
- `notebooks/` - Jupyter notebooks for exploratory analysis and figure generation (see `notebooks/README.md` for an overview)
- `resources/` - Gene sets, metadata, and supporting files


Each pipeline directory contains its own README with usage instructions.