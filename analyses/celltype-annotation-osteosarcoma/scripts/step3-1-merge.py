import scanpy as sc, anndata as ad, scrublet as scr
import scanpy.external as sce
import pandas as pd
from glob import glob
import numpy as np
from tqdm import tqdm
import os
from scipy.sparse import csr_matrix
import seaborn as sns
import matplotlib.pyplot as plt
import urllib.request
import decoupler as dc

import warnings
warnings.filterwarnings('ignore')

############################################################################################################

import importlib  
# import the previous step
step2_modules = importlib.import_module("step2-celltype-aucell")

def merge_and_process(
        matrix_folder = 'analyses/celltype-annotation-osteosarcoma/results/step1-processing',
        cell_type_folder = 'analyses/celltype-annotation-osteosarcoma/results/step2-processing',
        save_path = 'analyses/celltype-annotation-osteosarcoma/results/step3-processing'):

    #obtain all the adata objects
    all_adata = {}
    for f in tqdm(glob(matrix_folder + '/*.h5ad')):
        sample_id = f.split('/')[-1].replace('.h5ad', '')
        #read in the adata object
        adata = sc.read_h5ad(f)
        #read in cell type information and map onto the adata object
        cell_type = pd.read_csv(cell_type_folder + '/{s}.csv'.format(s = sample_id), index_col=0)
        cell_type['cell_type'] = cell_type.idxmax(axis = 1)
        cell_type_map = dict(zip(cell_type.index, cell_type['cell_type']))
        adata.obs['cell_type'] = adata.obs.index.map(cell_type_map)
        all_adata[sample_id] = adata
    all_adata = ad.concat(all_adata, label='sample', join = 'outer')
    all_adata.obs_names_make_unique()
    all_adata.var_names_make_unique()
    print('Cohort Size: ', all_adata.n_obs, all_adata.n_vars)

    #process the data
    all_adata = step2_modules.processing(all_adata)
    os.makedirs(save_path, exist_ok = True)
    #for some reason it doesn't like some of the column names that have `None` in them. we can drop them for now, just make a note of it.
    all_columns = all_adata.obs.columns
    all_adata.obs = all_adata.obs.dropna(axis=1, how='any')
    for col in all_adata.obs.columns:
        all_adata.obs[col] = all_adata.obs[col].astype(str)
    print('Dropped Columns: ', set(all_columns) - set(all_adata.obs.columns))
    all_adata.write_h5ad(save_path+'/all-samples.h5ad', compression = 'gzip')

    #save visualization
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(20, 10))
    sc.pl.umap(all_adata, color = ['sample_leiden_05'], ax = ax1)
    sc.pl.umap(all_adata, color = ['cell_type'], ax = ax2)
    plt.savefig(save_path + '/all-samples.png', dpi = 500)
    plt.close()

############################################################################################################

if __name__ == '__main__':
    merge_and_process()