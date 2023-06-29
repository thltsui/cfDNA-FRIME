# FRIME Simulation Repository

This repository contains scripts and data for generating simulations of FRIME (Fragmentation, Immigration, and Exit) processes, studying their stationary distributions, and comparing them with clinical cell-free DNA fragment profiles.

## Table of Contents

Description

Folder Structure

Usage

Data Preprocessing

Scripts

Simulation Results

Clinical Data Analysis Results


## Description

The aim of this repository is to provide a framework for generating simulations of FRIME processes and analyzing their stationary distributions. FRIME processes model the fragmentation, immigration, and exit of cell-free DNA fragments. The simulations are compared with clinical cell-free DNA fragment profiles to gain insights and understanding.

## Folder Structure

The repository is organized into the following folders:

**script**: Contains scripts related to the FRIME simulation and analysis.

**data_preprocessing**: Contains files and scripts for preprocessing the data used in the simulation and analysis.



After running the jupyter notebooks under the **script** folder, the following two folders will also be generated.

**simulation_results**: Stores the results and outputs generated from running the FRIME simulations.

**data_analysis**: Contains the analysis results and findings obtained from comparing the FRIME simulations with clinical cell-free DNA fragment profiles.

## Usage

To use this repository, follow these steps:

Clone the repository to your local machine:

```
git clone https://github.com/thltsui/cfDNA-FRIME.git

cd cfDNA-FRIME
```

Run the jupyter notebooks under **script**.


## Data Preprocessing

The Data Preprocessing folder contains the following files:


**script_1_data_preprocess.sh**: This script contains the codes and instructions to generate cell-free DNA fragment counts for our study.


**mtDNA_modeling_data.csv**: This CSV file is produced from **script_1_data_preprocess.sh** and contains cell-free DNA fragment counts for clinical data analysis.


## Scripts

The Scripts folder contains the following files:


**data_analysis.ipynb**: This Jupyter notebook is used to conduct cell-free DNA clinical data analysis.


**frime.py**: This Python script stores functions necessary to run FRIME simulations.


**simulation.ipynb**: This Jupyter notebook is used to run different numerical experiments on the stationary profile of FRIME simulations.


## Simulation Results

The **simulation_results** folder stores the results and outputs generated from running the FRIME simulations. The contents of this folder can be generated after running **simulation.ipynb**. In particular, the stationary profile plots stored in **simulation_results/numerical_experiment/clean_output** are used to produce Figure 3,4,5 in the manunuscript and Fig S1, S2, S9 in the supplementary materials.


## Clinical Data Analysis Results

The **data_analysis** folder contains the analysis results and findings obtained from comparing the FRIME simulations with clinical cell-free DNA fragment profiles. The content of this folder can be generated after running **data_analysis.ipynb**. They include plots of clinical cfDNA fragment profile for Fig 2,6,7,8 in the manuscript and Fig S3, S5-S8 in the supplementary materials, parameters of best-fitting FRIME simulations for Table 3 in the manuscript, and the $p$-value evolution of FRIME simulations for Fig S4.
