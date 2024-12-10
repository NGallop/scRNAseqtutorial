import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sc.settings.figdir = "plots/"

# read in raw data
raw_data = sc.read_h5ad('data/Hufallopiantube.h5ad')

# create a copy of the expression data matrix X in raw data to avoid changing original X
# create the layers dictionary to store alternative representations of X where "counts" is the key
raw_data.layers["counts"] = raw_data.X.copy()


# visualise the data in the terminal
print("Visualising the Observations in the AnnData object:")
print(raw_data.obs)
print("\n\n")
print("Highest Expressed Genes plot available. Please close to continue...")
#sc.pl.highest_expr_genes(raw_data, n_top=5)
print("\n\n")

####### Clean Data #######

# Filter out cells which have less than 100 genes expressed.
sc.pp.filter_cells(raw_data, min_genes= 100)

# Filter out genes which are expressed in less than 3 cells.
sc.pp.filter_genes(raw_data, min_cells= 3)

# Calculate quality control metrics. These are added to column of adata.obs or in adata.var_names.
sc.pp.calculate_qc_metrics(raw_data, percent_top=[50/100/200/500], inplace=True, log1p=False, )

# Use the above calculated QC metrics and visualise in a violin plot
print("Quality Control Metrics available. Please close to continue...")
#sc.pl.violin(raw_data, ['n_genes_by_counts','total_counts'], jitter=0.4, multi_panel=True)
print("\n\n")

# Filter again based on this to remove n_gene_by_counts of <6000 and 'total_counts' <1500000
raw_data = raw_data[raw_data.obs.n_genes_by_counts <6000,:]
data = raw_data[raw_data.obs.total_counts <1500000,:]
# Visualise again
print("Filtered Quality Control Metrics available. Please close to continue...")
#sc.pl.violin(data, ['n_genes_by_counts','total_counts'], jitter=0.4, multi_panel=True)
print("\n\n")


###### Dimensionality Reduction and Visualisation ######

'''
Steps to plot the data in a UMAP:
    1. Initial pre-processing steps
        - Normalise the counts per cell so each cell has 10000 counts
        - log transform the data
    2. Identify and crop the data to only the highly variable genes, then scale the data
    3. Compute the PCA
    4. Compute the Nearest Neighbour graph
    5. Compute the UMAP
'''

# 1a. Normalise the data
sc.pp.normalize_total(data, target_sum=1e4)

# 1b. Log transform the data
sc.pp.log1p(data)

# 2. Crop to highly variable genes
sc.pp.highly_variable_genes(data, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(data)
data.raw = data
data = data [:, data.var.highly_variable]
sc.pp.scale(data, max_value=10)

# 3. PCA
sc.tl.pca(data, svd_solver='arpack')
sc.pl.pca_variance_ratio(data, log=True)
sc.pl.pca(data, color=["Disease_stage"])

# 4. Nearest Neighbour
sc.pp.neighbors(data, n_pcs=10)

# 5. UMAP
sc.tl.umap(data)
sc.pl.umap(data, color=["Disease_stage"])



###### Leiden Clustering ######
'''
Cluster the cells to identify similar groups of cells

    By examining gene expression, determine what type of cells are in each cluster
    Using a dictionary, create a new obs. column assigning a label of cell type to leiden clusters
    Test the top differentially expressed genes in each leiden cluster/cell type
    Plot a heatmap of top differentially expressed genes in each leiden cluster/cell type
'''

# 1. Perform leiden clustering, then visualise clusters on the umap
sc.tl.leiden(data, resolution=0.2)
sc.pl.umap(data, color=['leiden'])

# 2. Calculate the top differentially expressed genes in each leiden cluster
sc.tl.rank_genes_groups(data,'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(data, n_genes=25)

# 3. Export the ranked genes with scores to a csv file
Export_csv = False
if Export_csv:
    topmarkers=data.uns['rank_genes_groups']
    groups = topmarkers['names'].dtype.names
    celltypemarkergenes = pd.DataFrame({group +'_'+key[:1]:topmarkers[key][group]
                for group in groups for key in ['names','scores']})

    celltypemarkergenes.to_csv('data/celltypemarkergenes_export.csv', index=False)

# 4. Plot the top 20 differentially expressed genes in a heatmap
sc.pl.rank_genes_groups_heatmap(data, n_genes=20, swap_axes=True, show_gene_labels=False, vmin=-3, vmax=3 )

# 5. Visualise marker genes to determine the cell types present in each cluster
gene_list = ["PAX8", "SOX17", "OVGP1", "MUC16", 
 #"FOXJ1", "CCDC78", "PIFO", 
 #"PTPRC", "CD3D", "CD14", 
 #"COL1A1", "COL3A1", "DCN", 
 "KRT17", "WFDC2"]
print("Final violin plot to decide dictionary. Close to close program...")
sc.pl.violin(data, gene_list, groupby = 'leiden', save='_gene_exp.png')


# 6. Create a dictionary defining leiden clusters and their respective cell types
celltypedict = {
        '0' : 'Secretory Epithelial-1',
        '1' : 'Secretory Epithelial-2',
        '2' : 'Ciliated Epithelial',
        '3' : 'STIC-1',
        '4' : 'STIC-2',
        '5' : 'Immune',
        '6' : 'Fibroblast'
        }

data.obs['Celltype'] = data.obs['leiden'].map(celltypedict)

# 7. Plot a UMPA coloured by leiden cell type clusters
sc.tl.umap(data)
sc.pl.umap(data, color=["Celltype"])