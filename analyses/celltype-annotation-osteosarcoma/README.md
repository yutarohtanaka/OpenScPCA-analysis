# Osteosarcoma (SCPCP000017) snRNA-seq Data Cell Type Annotation

## Description

Here, we perform cell type annotation for the samples in the Osteosarcoma snRNA-seq (SCPCP000017) dataset (n = 27) provided on the Alex's Lemonade Stand Open single-cell Pediatric Cancer Atlas (scPCA) portal. 
We describe and provide the scripts and intermediate files necessary to perform this analysis in this folder.  

## Usage

We assume that the data has first been downloaded as described below, and the working directory is at `~/OpenScPCA-analysis`. We recommend users perform this analysis in a virtual environment to prevent any library conflicts. 

0. install all requirements. 
    - `pip3 install -r analyses/celltype-annotation-osteosarcoma/requirements.txt` 
1. first, we remove all of the cells annotated as doublets. the output files will be stored as `results/step1-processing/{sample-id}.h5ad` files. the statistics for this step will be saved in `results/step1-processing/processing_stats.csv`.
    - `python3 analyses/celltype-annotation-osteosarcoma/scripts/step1-preprocessing.py`
2. 

## Input files

We use the steps to download the snRNA-seq expression matrices and other relevant data in the AnnData format into this workspace with the workflow described [here](https://openscpca.readthedocs.io/en/latest/getting-started/accessing-resources/getting-access-to-data/). Additionally, we download the doublet detection results that have been run using scDblFinder, described [here](https://openscpca.readthedocs.io/en/latest/getting-started/accessing-resources/getting-access-to-data/#accessing-scpca-module-results). 
```
./download-data.py --format AnnData --projects SCPCP000017
./download-results.py --modules doublet-detection --projects SCPCP000017
```
This data is downloaded into `data/2024-08-22/SCPCP000017` after this command is performed. 
Additionally, we provide a manually compiled dataset of marker genes (`celltype-annotation-osteosarcoma/marker-genes.csv`) that are used to score and annotate the nuclei. 

## Output files


## Software requirements

We perform this annotation in a Python (3.9.2) environment. All other libraries used in this analysis are listed in the (`celltype-annotation-osteosarcoma/requirements.txt`) file.

## Computational resources

We are performing this on our lab's Google Cloud environment with a Virtual Machine (e2-standard-32) instance. 
A smaller machine may be able to perform this analysis, but we have not tested this. 
