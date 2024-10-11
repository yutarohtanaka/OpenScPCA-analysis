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

def processing(adata, hvg = 3000):
    '''
    perform the processing steps on the adata object required to perform clustering and scoring.
    '''
    
    #save the raw counts in the counts layer.
    adata.layers["counts"] = adata.X.copy()
    #log-normalize counts to median cell count, and exclude highly expressed genes in this calculation.
    sc.pp.normalize_total(adata, target_sum=None, exclude_highly_expressed=True)
    sc.pp.log1p(adata)
    #save the log-normalized counts in the raw layer
    adata.raw = adata
    #select highly variable genes using the seurat v3 implementation.
    #we use 3000 genes as the default number of highly variable genes.
    sc.pp.highly_variable_genes(adata, n_top_genes = int(hvg), flavor = 'seurat_v3', layer='counts') 

    #before running PCA, scale the data.
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack", use_highly_variable = True, n_comps = 20)

    #carry out leiden clustering
    sc.pp.neighbors(adata, key_added = 'sample', use_rep = 'X_pca', n_neighbors = 20, n_pcs = 20)
    sc.tl.leiden(adata, key_added="sample_leiden_05", resolution=0.5, neighbors_key = 'sample')
    sc.tl.paga(adata, groups = 'sample_leiden_05', neighbors_key = 'sample')
    sc.pl.paga(adata, plot=False)
    sc.tl.umap(adata, neighbors_key = 'sample', min_dist = 0.5, init_pos='paga')

    return adata


def run_aucell(adata, cell_type_marker_df):
    '''
    score the nuclei for the cell types using the aucell method.
    '''
    adata.X = adata.layers['counts']
    dc.run_aucell(mat=adata, net = cell_type_marker_df, source = 'specific_cell_type', 
                  target = 'gene_id', min_n = 1, use_raw = False)
    df = adata.obsm['aucell_estimate']

    return df

def umap_visualization(adata, df, save_path):
    
    df['cell_type'] = df.idxmax(axis=1)
    adata.obs['cell_type'] = adata.obs.merge(df[['cell_type']], 
                                                    left_index = True, right_index = True)['cell_type']

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(15, 5))
    #ax1 = plt.subplot(111, frame_on=False)
    sc.pl.umap(adata, color = ['sample_leiden_05'], wspace = 0.25, ax = ax1)
    #ax1 = plt.subplot(111, frame_on=False)
    sc.pl.umap(adata, color = ['cell_type'], wspace = 0.25, ax = ax2)
    plt.savefig(save_path, dpi = 500)
    plt.close()

############################################################################################################

def scoring_cell_types(input_folder = 'analyses/celltype-annotation-osteosarcoma/results/step1-processing',
                       cell_type_marker_file = 'analyses/celltype-annotation-osteosarcoma/marker-genes.csv'):
    '''
    input is the folder path for the doublet-removed anndata object files.
    we've hard coded the path to the files for now, so make sure to change the path if necessary.
    '''
    #make sure you're in the right directory.
    assert os.getcwd().split('/')[-1] == 'OpenScPCA-analysis'

    #make the folder for the results
    out_folder = 'analyses/celltype-annotation-osteosarcoma/results/step2-processing'
    os.makedirs(out_folder, exist_ok=True)

    #load the cell type marker gene table.
    cell_type_marker_df = pd.read_csv(cell_type_marker_file)

    for adata_file in tqdm(glob(input_folder+'/*.h5ad')):
        sample_id = adata_file.split('/')[-1].replace('.h5ad', '')
        adata = sc.read_h5ad(adata_file)
        
        #run dimensionality reduction, and clustering on the data.
        adata = processing(adata)
        #run aucell to score cell types.
        celltype_scores = run_aucell(adata, cell_type_marker_df)
        celltype_scores.to_csv(out_folder + '/{s}.csv'.format(s = sample_id), index=True)

        #visualize the cell types on the UMAP.
        umap_visualization(adata, celltype_scores, save_path = out_folder + '/{s}.png'.format(s = sample_id))


############################################################################################################

if __name__ == '__main__':
    scoring_cell_types()