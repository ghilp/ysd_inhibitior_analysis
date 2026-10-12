# ysd_inhibitior_analysis
Analysis scripts for NGS sequencing data and AlphaFold 3 structural predictions supporting a yeast surface display platform for protease inhibitor screening.

## Overview

This repository contains Python scripts used to analyze sequencing data from selected yeast populations and evaluate structural confidence metrics for AlphaFold 3-predicted protein configurations.

The scripts support the analysis of protease inhibitor variants and the characterization of predicted protein complexes used in the study.

## Scripts
ngs_analysis.py

Processes FASTQ sequencing data and identifies amino acid sequences matching specified motifs or flanking sequences. The analysis can be used to identify and quantify sequence variants in selected populations.

alphafold3_plddt_analysis.py

Analyzes per-residue predicted Local Distance Difference Test (pLDDT) scores from AlphaFold 3 output files to summarize structural confidence over specified protein residue ranges.

alphafold3_pae_analysis.py

Analyzes predicted aligned error (PAE) values from AlphaFold 3 output files to characterize predicted confidence in the relative positioning of protein components.

## Requirements

Python 3 and the packages listed in requirements.txt are required.

Install the dependencies using:

pip install -r requirements.txt

## Input files

The scripts require input files appropriate to each analysis:

NGS analysis: FASTQ sequencing files, optionally compressed as .fastq.gz, depending on the input handling implemented in the script.
AlphaFold 3 analyses: Corresponding model .json and .cif files, as required by each script.

Input files are not included in this repository unless explicitly provided. Users must supply the appropriate files and configure the input paths and analysis parameters before running the scripts.

## Usage
1. Clone or download this repository.
2. Install the required Python packages.
3. Place the required input files in the desired directories.
4. Update the input paths and analysis parameters in the relevant script.
5. Run the scripts using Python.

The scripts may require configuration changes to match the user's input directory structure and analysis settings.

## Outputs

The scripts generate analysis outputs such as processed sequencing data, sequence-matching summaries, and structural confidence summaries. Exact output files and locations depend on the script configuration.

## Reproducibility

To reproduce the analyses reported in the associated study, users should use the appropriate input data and the same analysis parameters described in the manuscript and Supplementary Information.

## License

This project is distributed under the MIT License. See LICENSE for details.
