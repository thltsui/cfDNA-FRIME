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

**simulation_results**: Stores the results and outputs generated from running the FRIME simulations.

**data_analysis**: Contains the analysis results and findings obtained from comparing the FRIME simulations with clinical cell-free DNA fragment profiles.

## Usage

To use this repository, follow these steps:

Clone the repository to your local machine:

git clone https://github.com/thltsui/cfDNA-FRIME.git

cd cfDNA-FRIME

Run the jupyter notebooks under **script**.


## Data Preprocessing

The Data Preprocessing folder contains the following files:


**script_1_data_preprocess.sh**: This script contains the codes and instructions to generate cell-free DNA fragment counts for our study.


**mtDNA_modeling_data.csv**: This CSV file is produced from **script_1_data_preprocess.sh** and contains cell-free DNA fragment counts for clinical data analysis.


## Scripts

The Scripts folder contains the following files:


cfDNA_Data_analysis.ipynb: This Jupyter notebook is used to conduct cell-free DNA clinical data analysis.


Frag_Ex_Imm_Func.py: This Python script stores functions necessary to run FRIME simulations.


FRIME_Simulation.ipynb: This Jupyter notebook is used to run different numerical experiments on the stationary profile of FRIME simulations.


## Simulation Results

The Simulation Results folder stores the results and outputs generated from running the FRIME simulations. The contents of this folder can be generated after running **FRIME_Simulation.ipynb**.


## Clinical Data Analysis Results

The Clinical Data Analysis Results folder contains the analysis results and findings obtained from comparing the FRIME simulations with clinical cell-free DNA fragment profiles. The content of this folder can be generated after runnin **cfDNA_Data_analysis.ipynb**.
