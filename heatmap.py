import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
import os

def read_files_into_character(filepath):
    # Assuming the Excel files are in a tabular format
    return pd.read_excel(filepath).iloc[:, 0].tolist()

def convert_name(data):
    # Implement any specific name conversion if needed
    return data

# Reading data
AcuK_targets = read_files_into_character("/mnt/c/Users/59798/Documents/AcuK_M_Summary/noval_plot/study_iron_genes/AcuK_target_genes_orders_ace_v2.xls")
AcuK_activated = read_files_into_character("/mnt/c/Users/59798/Documents/For_paper2/acum_activated_genes.xls")
AcuK_repressed = read_files_into_character("/mnt/c/Users/59798/Documents/For_paper2/acum_repressed_genes.xls")
full_list_files = pd.read_csv("/mnt/c/Users/59798/Documents/For_paper2/ddsFull_fpkm_values_gene_v2_AcuK_AcuM_datasets.xls", sep="\t", index_col=0)
gene_list_central_carbon = convert_name(read_files_into_character("gene_list_central_pathways_full.xls"))

def plot_heatmaps_RNA_seq_acuK_acuM(gene_list_for_plot_normal):
    plt.figure(figsize=(10, 4.5))
    matrix_values = np.zeros((len(gene_list_for_plot_normal), 6))
    
    for i, gene in enumerate(gene_list_for_plot_normal):
        if gene in full_list_files.index:
            matrix_values[i, :] = full_list_files.loc[gene, [27, 30, 33, 11, 12, 13]]
    
    matrix_values[matrix_values == 0] = 1e-10
    ratio_matrix = np.log2(np.hstack((matrix_values[:, :3] / matrix_values[:, :3], matrix_values[:, 3:] / matrix_values[:, :3])))
    ratio_matrix[ratio_matrix > 4] = 4
    ratio_matrix[ratio_matrix < -4] = -4
    
    return ratio_matrix

ratio_matrix_plot = plot_heatmaps_RNA_seq_acuK_acuM(gene_list_central_carbon)

# Plot heatmap for RNA-seq data
plt.figure(figsize=(10, 10))
sns.heatmap(ratio_matrix_plot, cmap=sns.color_palette("coolwarm", as_cmap=True), center=0, cbar_kws={'label': 'Log2 Ratio'})
plt.show()

def plot_heatmaps_bindings_AcuK(gene_list_for_plot_normal):
    matrix_values = np.zeros((len(gene_list_for_plot_normal), 2))
    for i, gene in enumerate(gene_list_for_plot_normal):
        if gene in AcuK_targets:
            matrix_values[i, :] = 1
    return matrix_values

ratio_matrix_plot_bindings = plot_heatmaps_bindings_AcuK(gene_list_central_carbon)

# Plot heatmap for bindings data
plt.figure(figsize=(10, 10))
sns.heatmap(ratio_matrix_plot_bindings, cmap=sns.color_palette(["black", "yellow"]), cbar=False)
plt.show()

# Scatter plot for log2FoldChange vs -log10(padj)
plt.scatter(full_list_files['log2FoldChange'], -np.log10(full_list_files['padj']), c='gray', s=0.3)
repressed_genes = read_files_into_character("/Users/dongliguo/Dropbox/For_paper/acum_repressed_genes.xls")
activated_genes = read_files_into_character("/Users/dongliguo/Dropbox/For_paper/acum_activated_genes.xls")

plt.scatter(full_list_files.loc[repressed_genes, 'log2FoldChange'], 
            -np.log10(full_list_files.loc[repressed_genes, 'padj']), 
            c='skyblue', s=0.4, label='Repressed Genes')
plt.scatter(full_list_files.loc[activated_genes, 'log2FoldChange'], 
            -np.log10(full_list_files.loc[activated_genes, 'padj']), 
            c='tomato', s=0.4, label='Activated Genes')

plt.axvline(x=0.8, color='black')
plt.axvline(x=-0.8, color='black')
plt.ylim(0, 200)
plt.xlim(-15, 15)
plt.legend()
plt.show()

# Barplot
bar_data = [1876, 2434]
bar_labels = ['Label1', 'Label2']  # Add appropriate labels
plt.bar(bar_labels, bar_data, color=['orange', 'lightblue'])
plt.ylim(0, 3000)
plt.show()

# Read GFF file and write output
gff0_genome_location = pd.read_csv("/mnt/c/Users/59798/Documents/reference_sequences/A_nidulans/annotation/a_nidulans_exon_CDS_files_data_names.bed", sep="\t", header=None)
AcuK_genes_v1 = [gene for gene in AcuK_targets if gene in gff0_genome_location[3].tolist()]
pd.Series(AcuK_genes_v1).to_csv("AcuK_genes_v1.xls", index=False, header=False)

# Functions for make_tiles_for_promoter and normalization_tiles_for_promoter would need to be defined or imported

# Assuming you have these functions defined elsewhere
# tiles_2 = make_tiles_for_promoter_ATG("AcuK_genes_v1.xls")
# AcuK_glucose = normalization_tiles_for_promoter("/mnt/c/Users/59798/Documents/BAM_BED_files/AcuKmyc_Glu_TAGCTT_s_8_sequence.fastq.bam.bed")
# AcuK_acetate = normalization_tiles_for_promoter("/mnt/c/Users/59798/Documents/BAM_BED_files/AcuKmyc_Ace_GGCTAC_s_8_sequence.fastq.bam.bed")
# AcuK_proline = normalization_tiles_for_promoter("/mnt/c/Users/59798/Documents/BAM_BED_files/AcuKmyc_Pro_CTTGTA_s_8_sequence.fastq.gz.bam.bed")
# AcuM_glucose = normalization_tiles_for_promoter("/mnt/c/Users/59798/Documents/BAM_BED_files/AcuMHA_Glu_ATCACG_CW527-574_lib.sam.bam.bed")
# AcuM_acetate = normalization_tiles_for_promoter("/mnt/c/Users/59798/Documents/BAM_BED_files/AcuMHA_Ace_CGATGT_CW527-574_lib.sam.bam.bed")
# WT_HA_AcuK = normalization_tiles_for_promoter("/mnt/c/Users/59798/Documents/Bam_files_all/BAM_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_MH11036_CreA_1.bam.bed")

# dev.new(width=12, height=6.6)
plt.figure(figsize=(12, 6.6))
plt.subplots_adjust(wspace=0.5)
colors = sns.color_palette("YlOrBr", 200)
# Plotting heatmaps for promoters, assuming plot_heatmap_B is defined
# plot_heatmap_B(AcuK_glucose, AcuK_glucose, 1500, 321)
# plot_heatmap_B(AcuK_acetate, AcuK_acetate, 1500, 321)
# plot_heatmap_B(AcuM_glucose, AcuM_glucose, 1500, 321)
# plot_heatmap_B(AcuM_acetate, AcuM_acetate, 1500, 321)
# plot_heatmap_B(WT_HA_AcuK, WT_HA_AcuK, 1500, 321)

# Assuming normalization_tiles_for_gene_body_normal and other related functions are defined
# RNA_seq_summits_glucose = normalization_tiles_for_gene_body_normal("/Volumes/LiguoDisk/ribo_seq/16h_hypha_1_S23_L006.bam_summits_extsize_macs2_summits.bed")
# RNA_seq_summits_glucose = normalization_tiles_for_gene_body_normal("/Volumes/LiguoDisk/ribo_seq/16h_hypha_1_S23_L006.bam_summits_extsize_macs2_peaks_narrowPeaks.bed")
# plot(colMeans(RNA_seq_summits_glucose))
# plot_heatmap_B(RNA_seq_summits_glucose, RNA_seq_summits_glucose, 1500, 321)
# plot_heatmap_B(RNA_seq_summits_glucose, RNA_seq_summits_glucose, 1500, 321)
