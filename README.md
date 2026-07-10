# BioPython course

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Biopython](https://img.shields.io/badge/Biopython-1.81-4B8BBE?logo=python&logoColor=white)](https://biopython.org/)
[![Seaborn](https://img.shields.io/badge/Seaborn-Data_Viz-3776AB?logo=pandas&logoColor=white)](https://seaborn.pydata.org/)
[![PyDESeq2](https://img.shields.io/badge/PyDESeq2-RNA--seq-8A2BE2)](https://pydeseq2.readthedocs.io/)
[![Database](https://img.shields.io/badge/SQL-Relational_DB-003B57?logo=sqlite&logoColor=white)]()

> A repository dedicated to BioPython's class, in which numerous bioinformatic pipelines were explored to gain practical ability over them.

## Overview

This repository serves as a practical compilation of computational biology techniques implemented in Python. It is designed to demonstrate proficiency in handling biological data across different levels of abstraction, from parsing raw DNA sequences to executing statistical models for transcriptomics and managing structured data.

## Core Topics & Modules

### 1. DNA Sequence Handling (`Biopython`)
* **Sequence parsing & I/O:** Reading and writing standard formats (FASTA, GenBank) using `Bio.SeqIO`.
* **String manipulation:** Transcribing, translating, and calculating basic genomic metrics (e.g., GC content).
* **Alignment parsing:** Handling and extracting insights from sequence alignments.

### 2. Object-Oriented Programming (OOP) in Python
* **Biological Data Modeling:** Designing custom Python classes to represent biological entities (e.g., Genes, Transcripts, Proteins).
* **Encapsulation & Inheritance:** Building reusable, modular code pipelines that abstract away the complexity of data processing.
* **Magic Methods:** Implementing Pythonic protocols (`__init__`) for intuitive object interaction.

### 3. Data Visualization (`Seaborn`)
* **Exploratory Data Analysis (EDA):** Visualizing biological distributions and statistical trends.
* **Publication-Ready Plots:** Generating highly customizable visual outputs, including density plots, boxplots, and specialized genomic charts (e.g., Volcano plots, Heatmaps).

### 4. Differential Expression Analysis (`pydeseq2`)
* **RNA-seq Pipeline:** Running a full differential expression analysis natively in Python, eliminating the need to bridge with R.
* **Statistical Modeling:** Dispersions estimation, log-fold change calculations, and Wald tests to identify significantly altered gene expression profiles.

### 5. Relational Databases
* **Data Persistence:** Introduction to structured biological data storage.
* **SQL Queries:** Designing schemas, linking tables via foreign keys (e.g., mapping genes to pathways or functional annotations), and executing relational queries to extract meaningful biological subsets.

## Getting Started

### Prerequisites
To run the scripts and notebooks in this repository, you will need a working Python environment. The standard scientific stack is required.

```bash
# Recommended: Create a virtual environment via Conda or venv
conda create -n biopython_env python=3.10
conda activate biopython_env

# Install the required dependencies
pip install biopython seaborn pandas pydeseq2 numpy
```

## Repository structure

```text
.
├── data
│   ├── ConteosNormaliz.csv
│   ├── count_matrix.tsv
│   └── exp_dif.csv
├── doc
│   ├── Classes_1_2.ipynb
│   └── Sesion1.ipynb
├── LICENSE
├── results
├── src
│   ├── intro_db.ipynb
│   ├── intro_DEA.ipynb
│   ├── intro_seaborn.ipynb
│   ├── SalazarMendez_Pablo_Ejercicio2.py
│   ├── SalazarMendez_Pablo_Ejercicio3_4.py
│   ├── SalazarMendez_Pablo_Ejercicio5.py
│   └── sesion_2.ipynb
└── tmp
```
