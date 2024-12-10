import scanpy as sc
import numpy as np
import pandas as pd
import requests
import matplotlib.pyplot as plt

sc.settings.figdir = "plots/data_investigation/"

# read in raw data
raw_data = sc.read_h5ad('data/Hufallopiantube.h5ad')

# create a copy of the expression data matrix X in raw data to avoid changing original X
# create the layers dictionary to store alternative representations of X where "counts" is the key
raw_data.layers["counts"] = raw_data.X.copy()

print("Plot display suppressed; plots saved directly to file.")
print("\n\n")

# visualise the data in the terminal
print("Visualising the Observations in the AnnData object:")
print(raw_data.obs)
print("\n\n")
print("Highest Expressed Genes plot available.")
sc.pl.highest_expr_genes(raw_data, n_top=10, save='_highest_exp_genes.png', show=False)
print("\n\n")

####### Clean Data #######

# Filter out cells which have less than 100 genes expressed.
sc.pp.filter_cells(raw_data, min_genes= 100)

# Filter out genes which are expressed in less than 3 cells.
sc.pp.filter_genes(raw_data, min_cells= 3)

# Calculate quality control metrics. These are added to column of adata.obs or in adata.var_names.
sc.pp.calculate_qc_metrics(raw_data, percent_top=[50/100/200/500], inplace=True, log1p=False, )

# Use the above calculated QC metrics and visualise in a violin plot
print("Quality Control Metrics available.")
sc.pl.violin(raw_data, ['n_genes_by_counts','total_counts'], jitter=0.4, multi_panel=True, save='_raw_qc.png', show=False)
print("\n\n")

# Filter again based on this to remove n_gene_by_counts of <6000 and 'total_counts' <1500000
raw_data = raw_data[raw_data.obs.n_genes_by_counts <6000,:]
data = raw_data[raw_data.obs.total_counts <1500000,:]
# Visualise again
print("Filtered Quality Control Metrics available.")
sc.pl.violin(data, ['n_genes_by_counts','total_counts'], jitter=0.4, multi_panel=True, save='_sliced_qc.png', show=False)
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
sc.pl.highly_variable_genes(data, save='_high_variable_genes.png', show=False)
data.raw = data
data = data [:, data.var.highly_variable]
sc.pp.scale(data, max_value=10)

# 3. PCA
sc.tl.pca(data, svd_solver='arpack')
sc.pl.pca_variance_ratio(data, log=True, save="_variance_ratio.png", show=False)
sc.pl.pca(data, color=["Disease_stage"], save='_diseas_stage.png', show=False)

# 4. Nearest Neighbour
sc.pp.neighbors(data, n_pcs=8)

# 5. UMAP
sc.tl.umap(data)
sc.pl.umap(data, color=["Disease_stage"], save="_umap_disease_stage.png", show=False)



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
sc.pl.umap(data, color=['leiden'], save="_leiden_umap.png", show=False)

# 2. Calculate the top differentially expressed genes in each leiden cluster
sc.tl.rank_genes_groups(data,'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(data, n_genes=25, save="_ranked_genes.png", show=False)

# 3. Export the ranked genes with scores to a csv file
Export_csv = False
if Export_csv:
    topmarkers=data.uns['rank_genes_groups']
    groups = topmarkers['names'].dtype.names
    celltypemarkergenes = pd.DataFrame({group +'_'+key[:1]:topmarkers[key][group]
                for group in groups for key in ['names','scores']})

    celltypemarkergenes.to_csv('data/celltypemarkergenes_export.csv', index=False)

# 4. Plot the top 20 differentially expressed genes in a heatmap
sc.pl.rank_genes_groups_heatmap(data, n_genes=20, swap_axes=True, show_gene_labels=False, vmin=-3, vmax=3, save="_ranked_genes_groups.png", show=False)

# 5. Visualise marker genes to determine the cell types present in each cluster
sec_epi_list = ["PAX8", "SOX17", "OVGP1", "MUC16"]
cil_epi_list = ["FOXJ1", "CCDC78", "PIFO"]
immune_list = ["PTPRC", "CD3D", "CD14"]
fibroblast_list = ["COL1A1", "COL3A1", "DCN"]
stic_list = ["KRT17", "WFDC2"]

gene_lists = [sec_epi_list, cil_epi_list, immune_list, fibroblast_list, stic_list]

print("Final violin plot to decide dictionary. Close to close program...")
iteration = 0
for gene_list in gene_lists:
    iteration += 1
    if gene_list[0] == "PAX8":
        file_namer = "secretory_epithelial"
    elif gene_list[0] == "FOXJ1":
        file_namer = "ciliated_epithelial"
    elif gene_list[0] == "PTPRC":
        file_namer = "immune_cell"
    elif gene_list[0] == "COL1A1":
        file_namer = "fibroblast"
    elif gene_list[0] == "KRT17":
        file_namer = "STIC"
    else:
        file_namer = f"unknown{iteration}"
    sc.pl.violin(data, gene_list, groupby = 'leiden', save=f'_gene_exp_{file_namer}.png', show=False)


# 6. Create a dictionary defining leiden clusters and their respective cell types
celltypedict = {
        '0' : 'Secretory Epithelial-1',
        '1' : 'Secretory Epithelial-2',
        '2' : 'Ciliated Epithelial',
        '3' : 'STIC-1',
        '4' : 'Immune',
        '5' : 'STIC-2',
        '6' : 'Fibroblast'
        }

data.obs['Celltype'] = data.obs['leiden'].map(celltypedict)

# 7. Plot a UMPA coloured by leiden cell type clusters
sc.tl.umap(data)
sc.pl.umap(data, color=["Celltype"], save="_cell_type_umap.png", show=False)

###### GWAS associations ######

print("\n\n")
print("Requesting JSON from GWAS...")
response=requests.get("https://www.ebi.ac.uk/gwas/rest/api/studies/GCST003588/associations")

associations_data = response.json()['_embedded']['associations']
gwas = pd.json_normalize(associations_data, 'loci')
risk_alleles = gwas.explode('strongestRiskAlleles').explode('authorReportedGenes').reset_index(drop=True)

risk_alleles_df = pd.json_normalize(risk_alleles['strongestRiskAlleles'])
reportedgenes_df = pd.json_normalize(risk_alleles['authorReportedGenes'])

gwasgenes=pd.concat([risk_alleles_df, reportedgenes_df], axis=1)
gwasgenes = gwasgenes[['riskAlleleName', 'geneName']]

genes=[gene for gene in gwasgenes['geneName'].tolist() if gene !='Intergenic']

print("\n\n")
print("Dot plot of gene expression of carcinoma GWAS genes per leiden cluster")

missing_genes = [gene for gene in genes if gene not in data.raw.var_names]
valid_genes = [gene for gene in genes if gene in data.raw.var_names]

if missing_genes:
    print(f"WARNING: Missing genes: {missing_genes}. Removing form workflow.")  # Debug: list of genes not found
    sc.pl.dotplot(data, valid_genes, groupby='leiden', swap_axes=True, dendrogram=True, use_raw=True, save="_gwas.png", show=False)
else:
    sc.pl.dotplot(data, genes, groupby='leiden', swap_axes=True, dendrogram=True, use_raw=True, save="_gwas.png", show=False)

print("\n\n")
print("UMAP gene expression of the highly expressed carcinoa GWAS genes compared to the disease stage umap")
sc.pl.umap(data, color=['leiden','Disease_stage','TUBA1C','KRT8','HLA-C'], save="_gwas.png", show=False)