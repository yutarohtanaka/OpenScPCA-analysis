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

import warnings
warnings.filterwarnings('ignore')

############################################################################################################

def filter_doublets(adata_file, doublet_file, doublet_removed_adata):
    '''
    input is the file path for the raw data. we should open the original file, 
    '''
    #save stats of the doublet removal process.
    stats = [adata_file.split('/')[-1].split('_')[0]]

    adata = ad.read_h5ad(adata_file)
    doublets = pd.read_csv(doublet_file, sep = '\t').set_index('barcodes').rename(columns = {'score':'doublet_score', 'class':'doublet_class'})

    stats.extend([adata.n_obs, adata.n_vars]) #original file size
    adata.obs = adata.obs.merge(doublets, left_index=True, right_index=True, how='left')
    adata = adata[adata.obs['doublet_class'] == 'singlet']
    #in addition, remove all genes that are expressed in no cells after removal of doublets.
    sc.pp.filter_genes(adata, min_cells=1)
    stats.extend([adata.n_obs, adata.n_vars]) #file size post doublet removal

    adata.write(doublet_removed_adata, compression='gzip')

    return stats



############################################################################################################

def preprocessing(folder = 'data/2024-08-22/SCPCP000017'):
    '''
    input is the folder path for the raw data (eg. `data/2024-08-22/SCPCP000017`).
    we've hard coded the path to the files for now, so make sure to change the path if necessary.
    '''
    #make sure you're in the right directory.
    assert os.getcwd().split('/')[-1] == 'OpenScPCA-analysis'

    #make the folder for the results
    out_folder = 'analyses/celltype-annotation-osteosarcoma/results/step1-processing'
    os.makedirs(out_folder, exist_ok=True)

    stats_all = []

    for adata_file in tqdm(glob(folder+'/*/*_processed_rna.h5ad')):
        
        # format the path to the same sample's doublet detection file.
        doublet_file = adata_file.replace('_processed_rna.h5ad', '_processed_scdblfinder.tsv').replace('data/2024-08-22', 'data/2024-08-22/results/doublet-detection')
        # define output path for the doublet removed file.
        sample_id = adata_file.split('/')[-1].split('_')[0]
        doublet_removed_adata = out_folder + f'/{sample_id}.h5ad'

        # run the doublet removal function, and output stats. 
        stats = filter_doublets(adata_file, doublet_file, doublet_removed_adata)
        stats_all.append(stats)
    
    stats_all = pd.DataFrame(stats_all, columns = ['sample_id', 'n_obs_original', 'n_vars_original', 'n_obs_post_doublet_removal', 'n_vars_post_doublet_removal'])
    stats_all.to_csv(out_folder + '/processing_stats.csv', index=False)

if __name__ == '__main__':
    preprocessing()