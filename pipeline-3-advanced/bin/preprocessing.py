#!/usr/bin/python

"""
Preprocessing script for the scRNA-seq data.
"""

import argparse
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt


##

# Parse arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description='Preprocessing script for scRNA-seq data')
    
    # Required arguments
    parser.add_argument('--input', '-i', required=True,
                        help='Input directory path containing QC.h5ad file')
    
    # Optional arguments with defaults
    parser.add_argument('--n_hvgs', type=int, default=2000,
                        help='Number of highly variable genes (default: 2000)')
    parser.add_argument('--n_bins', type=int, default=20,
                        help='Number of bins for highly variable genes calculation (default: 20)')
    parser.add_argument('--n_pcs', type=int, default=30,
                        help='Number of principal components (default: 30)')
    parser.add_argument('--n_neighbors', type=int, default=15,
                        help='Number of neighbors for neighborhood graph (default: 15)')
    parser.add_argument('--metric', default='euclidean',
                        choices=['euclidean', 'manhattan', 'cosine', 'correlation'],
                        help='Distance metric for neighbor calculation (default: euclidean)')
    parser.add_argument('--target_sum', type=float, default=1e4,
                        help='Target sum for normalization (default: 1e4)')
    parser.add_argument('--min_umis', type=int, default=500,
                        help='Minimum number of UMIs per cell (default: 500)')
    parser.add_argument('--max_umis', type=int, default=np.inf,
                        help='Maximum number of UMIs per cell (default: 50000)')
    parser.add_argument('--min_genes', type=int, default=250,
                        help='Minimum number of genes per cell (default: 250)')
    parser.add_argument('--max_genes', type=int, default=np.inf,
                        help='Maximum number of genes per cell (default: 7500)')
    parser.add_argument('--max_percent_mt', type=float, default=0.15,
                        help='Maximum mitochondrial gene percentage (default: 20.0)')
    parser.add_argument('--gene_threshold', type=float, default=0.01,
                        help='Threshold for gene filtering as fraction of cells (default: 0.001)')
    parser.add_argument('--resolution', type=float, default=0.5,
                        help='Resolution parameter for Leiden clustering (default: 0.5)')
    
    return parser.parse_args()


##


# Parse arguments
args = parse_arguments()

# Assign parsed arguments to variables for compatibility
path_input = args.input
n_hvgs = args.n_hvgs
n_bins = args.n_bins
n_pcs = args.n_pcs
n_neighbors_1 = args.n_neighbors
metric = args.metric
target_sum = args.target_sum
min_umis = args.min_umis
max_umis = args.max_umis
min_genes = args.min_genes
max_genes = args.max_genes
max_percent_mt = args.max_percent_mt / 100.0  # Convert percentage to fraction
thr = args.gene_threshold
resolution = args.resolution

# Read data
adata = sc.read(path_input)

# Compute QC metrics
adata.obs['n_UMIs'] = adata.X.sum(axis=1).A1
adata.obs['n_genes'] = (adata.X > 0).sum(axis=1).A1
mt_counts = adata[:, adata.var_names.str.startswith('MT-')].X.sum(axis=1).A1
total_counts = adata.X.sum(axis=1).A1
adata.obs['percent_mt'] = mt_counts / total_counts

# Filter cells and genes
test = (
    (adata.obs['n_UMIs'] >= min_umis) & 
    (adata.obs['n_UMIs'] <= max_umis) & 
    (adata.obs['n_genes'] >= min_genes) & 
    (adata.obs['n_genes'] <= max_genes) & 
    (adata.obs['percent_mt'] <= max_percent_mt)
)
adata = adata[test, :].copy()

# Refilter genes, at least expressed in 0.1% of cells
cell_threshold = thr * adata.n_obs
test = (adata.X > 0).sum(axis=0).A1 >= cell_threshold
adata = adata[:, test].copy()

# Preprocess
sc.pp.normalize_total(adata, target_sum=target_sum)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=n_hvgs, n_bins=n_bins)
sc.tl.pca(adata, n_comps=n_pcs)
sc.pp.neighbors(adata, n_neighbors=n_neighbors_1, n_pcs=n_pcs, metric=metric)
sc.tl.umap(adata)

# Visualization
fig, ax = plt.subplots(1,3,figsize=(8,2.5))
sc.pl.pca(adata, color='n_UMIs', ax=ax[0], show=False)
sc.pl.pca(adata, color='n_genes', ax=ax[1], show=False)
sc.pl.pca(adata, color='percent_mt', ax=ax[2], show=False)
fig.tight_layout()
fig.savefig('QC_PCA.pdf', dpi=1000)

fig, ax = plt.subplots(1,3,figsize=(8,2))
sc.pl.umap(adata, color='n_UMIs', ax=ax[0], show=False)
sc.pl.umap(adata, color='n_genes', ax=ax[1], show=False)
sc.pl.umap(adata, color='sample', ax=ax[2], show=False, palette='Spectral')
fig.tight_layout()
fig.savefig('QC_UMAP.pdf', dpi=1000)


##


# Clustering
sc.tl.leiden(adata, resolution=resolution)

# Viz
fig, ax = plt.subplots(1,3,figsize=(8,2))
sc.pl.umap(adata, color='leiden', ax=ax[0], show=False)
fig.tight_layout()
fig.savefig('leiden.pdf', dpi=1000)