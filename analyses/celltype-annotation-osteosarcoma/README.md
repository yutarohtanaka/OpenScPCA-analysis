# Osteosarcoma (SCPCP000017) snRNA-seq Data Cell Type Annotation

## Description

Here, we perform cell type annotation for the samples in the Osteosarcoma snRNA-seq (SCPCP000017) dataset (n = 27) provided on the Alex's Lemonade Stand Open single-cell Pediatric Cancer Atlas (scPCA) portal. 
We describe and provide the scripts and intermediate files necessary to perform this analysis in this folder. We've aimed to make the scripts as transferable to other datasets as possible. 

## Usage

We assume that the data has first been downloaded as described below, and the working directory is at `~/OpenScPCA-analysis`. We recommend users perform this analysis in a virtual environment to prevent any library conflicts. 

0. **Environment Prep**
install all requirements. 
    - **`pip3 install -r analyses/celltype-annotation-osteosarcoma/requirements.txt`**
1. **Doublet Removal**
first, we remove all of the cells annotated as doublets. the output files will be stored as `results/step1-processing/{sample-id}.h5ad` files. the statistics for this step will be saved in `results/step1-processing/processing_stats.csv`.
    - **`python3 analyses/celltype-annotation-osteosarcoma/scripts/step1-preprocessing.py`**
2. **Scoring Cell Type Gene Sets**
next, we run AUCell per sample on the cell type marker gene sets to score the nuclei on the cell types we broadly expect to see in this dataset. we use the `results/step1-processing/{sample-id}.h5ad` files as input, and we will score these and save them as `results/step2-processing/{sample-id}.csv` files. additionally, we generate a `results/step2-processing/{sample-id}.png` file with the UMAP visualization of the sample, coloured with clustering performed and scoring-based cell types.  
    - **`python3 analyses/celltype-annotation-osteosarcoma/scripts/step2-celltype-aucell.py`**
3. **Orthogonal Validation of Cell Type Annotations**

4. **Annotation of Malignant Cells**
5. **Post-Processing, Compilation**

## Input files

We use the steps to download the snRNA-seq expression matrices and other relevant data in the AnnData format into this workspace with the workflow described [here](https://openscpca.readthedocs.io/en/latest/getting-started/accessing-resources/getting-access-to-data/). Additionally, we download the doublet detection results that have been run using scDblFinder, described [here](https://openscpca.readthedocs.io/en/latest/getting-started/accessing-resources/getting-access-to-data/#accessing-scpca-module-results). 
```
./download-data.py --format AnnData --projects SCPCP000017
./download-results.py --modules doublet-detection --projects SCPCP000017
```
This data is downloaded into `data/2024-08-22/SCPCP000017` after this command is performed. 
Additionally, we provide a manually compiled dataset of marker genes (`celltype-annotation-osteosarcoma/marker-genes.csv`) that are used to score and annotate the nuclei. 

## Output files

We have intentionally made incremental saves after each step of the cell type annotation analysis workflow described here, to make it easier to rerun from intermediate steps if necessary.

```
results
├── step1-processing
│   ├── {sample-id}.h5ad #contains the raw AnnData object with detected doublets removed.
│   └── processing_stats.csv #pre/post doublet removal stats.
├── step2-processing
│   ├── {sample-id}.csv #contains the cell type scores and inferred cell type.  
│   ├── {sample-id}.png #UMAP visualization w/ clustering and scoring-based cell type annotations.
└── README.md
```

## Software requirements

We perform this annotation in a Python (3.9.2) environment. All other libraries used in this analysis are listed in the (`celltype-annotation-osteosarcoma/requirements.txt`) file.

## Computational resources

We are performing this on our lab's Google Cloud environment with a Virtual Machine (e2-standard-32) instance. 
A smaller machine may be able to perform this analysis, but we have not tested this. 
