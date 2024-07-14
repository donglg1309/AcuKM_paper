# Load necessary functions
source("/mnt/c/Users/59798/Documents/ANAlyses_folder/AcuK_M_FacB_Summary/FacB_analysis/make_tiles_functions.R")

# Read data
AcuK_targets <- read_files_into_character("/mnt/c/Users/59798/Documents/AcuK_M_Summary/noval_plot/study_iron_genes/AcuK_target_genes_orders_ace_v2.xls")
AcuK_activated <- read_files_into_character("/mnt/c/Users/59798/Documents/For_paper2/acum_activated_genes.xls")
AcuK_repressed <- read_files_into_character("/mnt/c/Users/59798/Documents/For_paper2/acum_repressed_genes.xls")
full_list_files <- read.table("/mnt/c/Users/59798/Documents/For_paper2/ddsFull_fpkm_values_gene_v2_AcuK_AcuM_datasets.xls", sep="\t", row.names=1, header=TRUE)
gene_list_central_carbon <- comvert_name(read_files_into_character("gene_list_central_pathways_full.xls"))

# Function to plot heatmaps for RNA-seq data
plot_heatmaps_RNA_seq_acuK_acuM <- function(gene_list_for_plot_normal) {
  dev.new(width=10, height=4.5)
  matrix_values <- matrix(0, nrow=length(gene_list_for_plot_normal), ncol=6)
  
  for (i in 1:length(gene_list_for_plot_normal)) {
    matrix_values[i, ] <- as.numeric(full_list_files[rownames(full_list_files) %in% gene_list_for_plot_normal[i], c(27, 30, 33, 11, 12, 13)])
  }
  
  matrix_values[matrix_values == 0] <- 0.0000000001
  ratio_matrix <- log2(cbind(matrix_values[, 1:3] / matrix_values[, 1:3], matrix_values[, 4:6] / matrix_values[, 1:3]))
  ratio_matrix[ratio_matrix > 4] <- 4
  ratio_matrix[ratio_matrix < -4] <- -4
  
  return(as.matrix(ratio_matrix))
}

# Plot heatmap for RNA-seq data
my_palette <- colorRampPalette(c("blue", "blue", "blue", "whitesmoke", "whitesmoke", "red", "red", "red"))(n = 100)
ratio_matrix_plot <- plot_heatmaps_RNA_seq_acuK_acuM(gene_list_central_carbon)
heatmap.2(ratio_matrix_plot, symm=T, scale='none', symkey=T, symbreaks=T, sepcolor="black", trace="none", cexRow=1.3, density.info="none", col=my_palette, Colv=FALSE, Rowv=FALSE, breaks=seq(-4, 4, length.out=101), colsep=0:ncol(as.matrix(ratio_matrix_plot)), rowsep=0:nrow(as.matrix(ratio_matrix_plot)), sepwidth=c(0.0001, 0.00001))

# Function to plot heatmaps for AcuK bindings
plot_heatmaps_bindings_AcuK <- function(gene_list_for_plot_normal) {
  matrix_values <- matrix(0, nrow=length(gene_list_for_plot_normal), ncol=2)
  matrix_values[gene_list_for_plot_normal %in% AcuK_targets, ] <- 1
  return(as.matrix(matrix_values))
}

# Plot heatmap for AcuK bindings
my_palette <- colorRampPalette(c("black", "black", "yellow", "yellow"))(n = 100)
ratio_matrix_plot <- plot_heatmaps_bindings_AcuK(gene_list_central_carbon)
heatmap.2(ratio_matrix_plot, symm=T, scale='none', symkey=T, symbreaks=T, sepcolor="black", trace="none", cexRow=1.3, density.info="none", col=my_palette, Colv=FALSE, Rowv=FALSE, breaks=seq(0, 1, length.out=101), colsep=c(0, 0.5), rowsep=0:nrow(as.matrix(ratio_matrix_plot)), sepwidth=c(0.0001, 0.0001))

# Plot log2FoldChange vs -log10(padj)
plot(full_list_files$log2FoldChange, -log10(full_list_files$padj), pch=15, cex=0.3, ylim=c(0, 200), xlim=c(-15, 15))
points(
  full_list_files$log2FoldChange[as.character(full_list_files$Row.names) %in% read_files_into_character("/Users/dongliguo/Dropbox/For_paper/acum_repressed_genes.xls")],
  -log10(full_list_files$padj)[as.character(full_list_files$Row.names) %in% read_files_into_character("/Users/dongliguo/Dropbox/For_paper/acum_repressed_genes.xls")],
  pch=15, col='skyblue2', lwd=1, cex=0.4
)
points(
  full_list_files$log2FoldChange[as.character(full_list_files$Row.names) %in% read_files_into_character("/Users/dongliguo/Dropbox/For_paper/acum_activated_genes.xls")],
  -log10(full_list_files$padj)[as.character(full_list_files$Row.names) %in% read_files_into_character("/Users/dongliguo/Dropbox/For_paper/acum_activated_genes.xls")],
  pch=15, col='tomato', lwd=1, cex=0.4
)
abline(v=0.8, col="black")
abline(v=-0.8, col="black")

# Barplot
barplot(c(1876, 2434), col=c("orange", "lightblue"), ylim=c(0, 3000))

# Read GFF genome location data
gff0_genome_location <- read.table("/mnt/c/Users/59798/Documents/reference_sequences/A_nidulans/annotation/a_nidulans_exon_CDS_files_data_names.bed")
write.table(AcuK_targets[AcuK_targets %in% as.character(gff0_genome_location$V4)], file="AcuK_genes_v1.xls", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Normalization tiles for promoter data
tiles_2 <- make_tiles_for_promoter_ATG("AcuK_genes_v1.xls")
AcuK_glucose <- normalization_tiles_for_promoter("/mnt/c/Users/59798/Documents/BAM_BED_files/AcuKmyc_Glu_TAGCTT_s_8_sequence.fastq.bam.bed")
AcuK_acetate <- normalization_tiles_for_promoter("/mnt/c/Users/59798/Documents/BAM_BED_files/AcuKmyc_Ace_GGCTAC_s_8_sequence.fastq.bam.bed")
AcuK_proline <- normalization_tiles_for_promoter("/mnt/c/Users/59798/Documents/BAM_BED_files/AcuKmyc_Pro_CTTGTA_s_8_sequence.fastq.gz.bam.bed")
AcuM_glucose <- normalization_tiles_for_promoter("/mnt/c/Users/59798/Documents/BAM_BED_files/AcuMHA_Glu_ATCACG_CW527-574_lib.sam.bam.bed")
AcuM_acetate <- normalization_tiles_for_promoter("/mnt/c/Users/59798/Documents/BAM_BED_files/AcuMHA_Ace_CGATGT_CW527-574_lib.sam.bam.bed")
AcuM_proline <- normalization_tiles_for_promoter("/mnt/c/Users/59798/Documents/BAM_BED_files/")
WT_HA_AcuK <- normalization_tiles_for_promoter("/mnt/c/Users/59798/Documents/Bam_files_all/BAM_files/CreA_new_ChIP_seqs_files/bam_files_out/GLU_MH11036_CreA_1.bam.bed")

# Plot heatmaps for various conditions
dev.new(width=12, height=6.6)
par(mar=c(1, 1, 1, 1), mfrow=c(1, 6), xpd=FALSE)
colors <- colorRampPalette(c('black', 'black', 'yellow', 'yellow', 'yellow3', 'yellow3', 'orange', 'orange'))(200)
plot_heatmap_B(AcuK_glucose, AcuK_glucose, 1500, 321)
plot_heatmap_B(AcuK_acetate, AcuK_acetate, 1500, 321)
plot_heatmap_B(AcuM_glucose, AcuM_glucose, 1500, 321)
plot_heatmap_B(AcuM_acetate, AcuM_acetate, 1500, 321)
plot_heatmap_B(WT_HA_AcuK, WT_HA_AcuK, 1500, 321)

# Normalization tiles for gene body
tiles_2 <- make_tiles_for_body_with_normal_promoters("/Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/a_nidulans_exon_CDS_files_data_names.xls")
RNA_seq_summits_glucose <- normalization_tiles_for_gene_body_normal("/Volumes/LiguoDisk/ribo_seq/16h_hypha_1_S23_L006.bam_summits_extsize_macs
