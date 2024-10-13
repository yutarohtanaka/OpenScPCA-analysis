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
import infercnvpy as cnv
import wget 
import warnings
warnings.filterwarnings('ignore')

############################################################################################################

def compile(
        normal_cell_annot_folder = 'analyses/celltype-annotation-osteosarcoma/results/step2-processing',
        infercnv_annot_file = 'analyses/celltype-annotation-osteosarcoma/results/step4-processing/all_samples_malignant_cell_annotation.csv',
        out_folder = 'analyses/celltype-annotation-osteosarcoma/results/celltype-annotations',
        normal_cell_types = ['B', 'Endothelial', 'NK/T', 'Epithelial', 'Myeloid', 'Osteoclast']):
    
    #make sure the output folder exists
    os.makedirs(out_folder, exist_ok = True)

    #read in the normal cell annotations
    df_all = pd.DataFrame()
    for f in glob(normal_cell_annot_folder+'/*.csv'):
        df = pd.read_csv(f)
        df['sample'] = os.path.basename(f).split('.')[0]
        df_all = pd.concat([df_all, df])
    df_all = df_all.set_index('Unnamed: 0')
    df_all['cell_type'] = df_all.drop(columns = ['sample']).idxmax(axis = 1)
    df_all = df_all[df_all.cell_type.isin(normal_cell_types)]
    df_all = df_all[['sample', 'cell_type']]

    #read in the infercnv annotations
    malignant_annot = pd.read_csv(infercnv_annot_file)
    malignant_annot['cell_type'] = malignant_annot.apply(lambda x: 'Tumor' if x['malignant'] == 'malignant' else x['cell_type'], 
                                                         axis = 1)
    malignant_annot = malignant_annot[['sample', 'cell_type']]
    
    #merge the two annotations
    df_all = pd.concat([df_all, malignant_annot])
    df_all.to_csv(out_folder + '/all-samples-annotations.csv', index=True)

############################################################################################################

if __name__ == '__main__':
    compile()
