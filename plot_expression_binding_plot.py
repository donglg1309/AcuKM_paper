import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Function to read files into a list of strings
def read_files_into_list(file_path):
    return pd.read_csv(file_path, header=None).iloc[:, 0].tolist()

# Function to create a matrix for RNA-seq heatmap
def create_rna_seq_matrix(gene_list, full_list_df):
    matrix_values = np.zeros((len(gene_list), 6))
    for i, gene in enumerate(gene_list):
        if gene in full_list_df.index:
            matrix_values[i, :] = full_list_df.loc[gene, [27, 30, 33, 11, 12, 13]].values
    matrix_values[matrix_values == 0] = 1e-10  # Avoid division by zero
    ratio_matrix = np.log2(np.c_[matrix_values[:, :3] / matrix_values[:, :3], matrix_values[:, 3:] / matrix_values[:, :3]])
    ratio_matrix = np.clip(ratio_matrix, -4, 4)
    return ratio_matrix

# Function to create a matrix for binding heatmap
def create_binding_matrix(gene_list, target_genes):
    matrix_values = np.zeros((len(gene_list), 2))
    matrix_values[[gene in target_genes for gene in gene_list], :] = 1
    return matrix_values

# Function to plot a heatmap
def plot_heatmap(matrix_plot, palette, vmin, vmax, file_name):
    plt.figure(figsize=(10, len(matrix_plot) * 0.5))
    sns.heatmap(matrix_plot, cmap=palette, cbar=True, vmin=vmin, vmax=vmax, linewidths=0.1, linecolor='black')
    plt.savefig(file_name)
    plt.close()

# Load the data
AcuK_targets = read_files_into_list("/mnt/c/Users/59798/Documents/AcuK_M_Summary/noval_plot/study_iron_genes/AcuK_target_genes_orders_ace_v2.xls")
AcuK_activated = read_files_into_list("/mnt/c/Users/59798/Documents/For_paper2/acum_activated_genes.xls")
AcuK_repressed = read_files_into_list("/mnt/c/Users/59798/Documents/For_paper2/acum_repressed_genes.xls")
full_list_df = pd.read_csv("/mnt/c/Users/59798/Documents/For_paper2/ddsFull_fpkm_values_gene_v2_AcuK_AcuM_datasets.xls", sep="\t", index_col=0)

# Central carbon metabolism analysis
gene_list_central_carbon = pd.read_csv("amino_acid_heatmap_plot.xls", sep=",")['gene_name'].unique()
pd.Series(gene_list_central_carbon).to_csv("gene_list_central_carbon_plot.xls", index=False, header=False)

# Plot RNA-seq heatmap for central carbon metabolism
ratio_matrix_plot = create_rna_seq_matrix(gene_list_central_carbon, full_list_df)
my_palette = sns.color_palette("coolwarm", as_cmap=True)
plot_heatmap(ratio_matrix_plot, my_palette, -4, 4, "central_carbon_rna_seq_heatmap.png")

# Plot binding heatmap for central carbon metabolism
ratio_matrix_plot = create_binding_matrix(gene_list_central_carbon, AcuK_targets)
my_palette = sns.color_palette("YlOrBr", as_cmap=True)
plot_heatmap(ratio_matrix_plot, my_palette, 0, 1, "central_carbon_binding_heatmap.png")

# Secondary metabolism analysis
gene_list_secondary_metabolism = pd.read_csv("secondary_metabolism_heatmap_plot.xls", sep=",")['gene_name'].unique()
pd.Series(gene_list_secondary_metabolism).to_csv("secondary_metabolism_heatmap_plot_v2.xls", index=False, header=False)

# Plot RNA-seq heatmap for secondary metabolism
ratio_matrix_plot = create_rna_seq_matrix(gene_list_secondary_metabolism, full_list_df)
my_palette = sns.color_palette("coolwarm", as_cmap=True)
plot_heatmap(ratio_matrix_plot, my_palette, -4, 4, "secondary_metabolism_rna_seq_heatmap.png")

# Plot binding heatmap for secondary metabolism
ratio_matrix_plot = create_binding_matrix(gene_list_secondary_metabolism, AcuK_targets)
my_palette = sns.color_palette("YlOrBr", as_cmap=True)
plot_heatmap(ratio_matrix_plot, my_palette, 0, 1, "secondary_metabolism_binding_heatmap.png")
