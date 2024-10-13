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

def map_chr_coordinates(adata, out_folder = 'analyses/celltype-annotation-osteosarcoma/results/step4-processing',):
    '''
    for all genes in the adata, map the chromosome coordinates to the hg38 genome. 
    '''
    #download the chromosome coordinates and process the file
    chr_coordinates_files = out_folder + '/Homo_sapiens.GRCh38.112.gtf.gz'
    if os.path.exists(chr_coordinates_files) == False:
        wget.download('https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz', chr_coordinates_files)
    gtf = pd.read_csv(chr_coordinates_files, delimiter = '\t', comment='#', header = None)
    gtf.columns = ['chromosome', "source", "feature", "start", "end", 'score', 'strand', 'phase', 'attributes']
    gtf = gtf[gtf['feature'] == 'gene']
    gtf['attributes'] = gtf.attributes.apply(lambda x: {i.split(' ')[0]:i.split(' ')[1].replace('"','') for i in x.strip(';').split('; ')})
    gtf['gene_id'] = gtf.attributes.apply(lambda x: x['gene_id'])
    gene_ids_adata = adata.var_names.tolist()
    gtf = gtf.loc[gtf.gene_id.isin(gene_ids_adata), :]
    gene_pos = gtf[['chromosome', 'start', 'end', 'gene_id']].drop_duplicates(subset='gene_id')
    
    #map it onto the anndata object
    adata = ad.AnnData(X = adata.layers['counts'], obs = adata.obs,
                       var = adata.var.merge(gene_pos, left_index = True, 
                                             right_on = 'gene_id', how = 'left').set_index('gene_id'))
    #drop any unmapped genes
    adata = adata[:, adata.var.chromosome.notna()]
    #drop any non-autosomal genes
    adata.var['chromosome'] = adata.var['chromosome'].astype(str)
    adata = adata[:, adata.var.chromosome.isin([str(i) for i in range(1,23)])]
    adata.var['chromosome'] = adata.var['chromosome'].apply(lambda x: 'chr' + x)

    return adata 

def downsample(adata, out_path, 
               cell_type_str = 'cell_type', normal_cell_types = [], fraction = 0.10, random_state = 42):
    '''
    create an object of nuclei annotated as normal cells. 
    in this cohort, we will downsample this object to 10% of all non-mesenchymal cells in the cohort.
    for reproducibility, save this object. 
    ''' 
    normal_cohort = adata[adata.obs[cell_type_str].isin(normal_cell_types)].copy()
    normal_cohort_subsample = sc.pp.subsample(normal_cohort, fraction = 0.10, random_state = 42, copy = True)
    normal_cohort_subsample.obs['subsampled'] = 'normal_subsample'
    normal_cohort_subsample.write_h5ad(out_path, compression = 'gzip')

    return normal_cohort_subsample

def infercnv(adata_sample, subsampled_normal_pop, out_path, normal_cell_types, cell_type_str = 'cell_type', ):
    '''
    for each sample, we take the mesenchymal cells, and add in the "normal" control cohort. 
    then, we run infercnvpy, and save the anndata objects per-sample. we've hard encoded the parameters here, but they may be tweaked to adjust for other cohorts.
    '''
    #subset the mesenchymal cells
    adata_sample = adata_sample[adata_sample.obs[cell_type_str].isin(normal_cell_types) == False].copy()
    adata_sample.obs['subsampled'] = 'mesenchymal_sample'
    #add in the normal control cohort
    adata_sample = adata_sample.concatenate(subsampled_normal_pop, index_unique = None)
    if 'chromosome-0' in adata_sample.var.columns:
        adata_sample.var['chromosome'] = adata_sample.var['chromosome-0']
    
    ref_cells = subsampled_normal_pop.obs[cell_type_str].unique().tolist()
    cnv.tl.infercnv(adata_sample, reference_key = cell_type_str, lfc_clip = 0.05, 
                    reference_cat = ref_cells, window_size = 250, step = 1)
    cnv.tl.pca(adata_sample, svd_solver="arpack")
    cnv.pp.neighbors(adata_sample, n_neighbors = 20, n_pcs = 20)
    cnv.tl.leiden(adata_sample, resolution=1)
    cnv.tl.umap(adata_sample, min_dist = 0.5)
    cnv.tl.cnv_score(adata_sample)

    adata_sample.write_h5ad(out_path, compression = 'gzip')

def determine_malignant_cells(adata_sample, normal_cell_types, cell_type_str = 'cell_type'):
    '''
    determine the malignant cells in the cohort. as osteosarcoma does not have a ubiquitous copy number alteration we can use as a marker, 
    we will use the mean + 1 standard deviation of the cnv_score in the normal control cohort as a cutoff.
    '''
    cutoff = adata_sample.obs[adata_sample.obs.subsampled == 'normal_subsample'].cnv_score.mean() + \
        adata_sample.obs[adata_sample.obs.subsampled == 'normal_subsample'].cnv_score.std()
    
    malignant_annot = adata_sample.obs[adata_sample.obs.subsampled == 'mesenchymal_sample'][['cnv_score', cell_type_str]]
    malignant_annot['malignant'] = malignant_annot['cnv_score'].apply(lambda x: 'malignant' if x > cutoff else 'normal')
    
    return malignant_annot
    
############################################################################################################

def run_infercnv(adata_path = 'analyses/celltype-annotation-osteosarcoma/results/step3-processing/all-samples.h5ad',
                 out_folder = 'analyses/celltype-annotation-osteosarcoma/results/step4-processing'):
    os.makedirs(out_folder, exist_ok = True)
    #read in anndata object
    adata = sc.read_h5ad(adata_path)
    #first, map all genes in the matrix to chromosome coordinates.
    adata = map_chr_coordinates(adata)

    #make the downsampled "control" population
    normal_cell_types = ['B', 'Endothelial', 'NK/T', 'Epithelial', 'Myeloid', 'Osteoclast']
    subsampled_normal_pop = downsample(adata, cell_type_str='cell_type', 
                                       normal_cell_types = normal_cell_types,
                                       out_path = out_folder + '/normal-subsampled-cells.h5ad')
    
    #run infercnv per-sample
    for sample in tqdm(adata.obs['sample'].unique()):
        adata_sample = adata[adata.obs['sample'] == sample].copy()
        out_path = out_folder + '/' + sample + '_infercnv_out.h5ad'
        infercnv(adata_sample, subsampled_normal_pop, out_path, normal_cell_types, cell_type_str = 'cell_type', )
    
    #determine the malignant cells in the cohort
    malignant_cell_annotation_df = pd.DataFrame()
    for f in tqdm(glob(out_folder + '/*infercnv_out.h5ad')):
        sample_id = f.split('/')[-1].replace('_infercnv_out.h5ad', '')
        adata = sc.read_h5ad(f)
        
        malignant_annot = determine_malignant_cells(adata, normal_cell_types, cell_type_str = 'cell_type')
        malignant_annot['sample'] = sample_id
        malignant_cell_annotation_df = pd.concat([malignant_cell_annotation_df, malignant_annot])
    
    malignant_cell_annotation_df.to_csv(out_folder + '/all_samples_malignant_cell_annotation.csv', index = True)
    
    print('We identify {n} malignant cells in the cohort.'.format(n = malignant_cell_annotation_df[malignant_cell_annotation_df.malignant == 'malignant'].shape[0]))

############################################################################################################

if __name__ == '__main__':
    run_infercnv()
