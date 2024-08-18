import numpy as np
import pandas as pd
import pybedtools
import pyBigWig
import matplotlib.pyplot as plt
from pybedtools import BedTool
from matplotlib.colors import LinearSegmentedColormap

def make_tiles_for_summit(bed_file, seq_range, chromatin_length_file):
    # Import BED file using pybedtools
    peaks = BedTool(bed_file)
    
    # Create tiles around peaks
    tiles = []
    for feature in peaks:
        center = (int(feature.end))
        tiles.extend([center + x for x in np.linspace(seq_range[0], seq_range[2], 301)])
    
    # Load chromatin length information
    chromatin_length = pd.read_table(chromatin_length_file, header=None)
    seqinfo = dict(zip(chromatin_length[0], chromatin_length[1]))
    
    # Create GRanges equivalent in Python
    tile_bed = []
    for i, start in enumerate(tiles):
        chrom = feature.chrom  # Assuming all peaks are from the same chromosome
        tile_bed.append([chrom, start, start + 50])
    
    # Convert to BedTool object
    tile_bedtool = BedTool(tile_bed)
    return tile_bedtool

def normalization_tiles_for_summits(bed_file, tiles, seq_range, genome_size_file):
    # Import BED file using pybedtools
    peaks = BedTool(bed_file)
    
    # Adjust start and width based on strand (dummy implementation as strands are not handled)
    adjusted_peaks = []
    for feature in peaks:
        if feature.strand == '-':
            feature.start -= (200 - (feature.end - feature.start))
        feature.end = feature.start + 200
        adjusted_peaks.append(feature)
    
    # Convert back to BedTool object
    adjusted_bed = BedTool(adjusted_peaks)
    
    # Count overlaps between tiles and peaks
    overlap_counts = tiles.intersect(adjusted_bed, c=True).to_dataframe()['score'].values
    
    # Reshape to matrix form
    num_peaks = len(overlap_counts) // 301
    overlap_matrix = overlap_counts.reshape(num_peaks, 301)
    
    # Normalize the matrix
    total_peaks = len(adjusted_bed)
    norm_matrix = overlap_matrix * 1000000 / (0.05 * total_peaks)
    
    return norm_matrix

def plot_heatmap(matrix, max_density, k):
    # Clip values above max_density
    matrix[matrix > max_density] = max_density
    
    # Plot heatmap
    colors = LinearSegmentedColormap.from_list("custom_cmap", ['black', 'yellow', 'orange'])
    plt.imshow(matrix, aspect='auto', cmap=colors, interpolation='none', vmin=0, vmax=max_density)
    plt.axis('off')
    plt.show()

def plot_figures(summit_file, chip_seq_file, max_density):
    # Define constants
    seq_range = [-1500, 1500]
    chromatin_length_file = "/mnt/c/Users/59798/Documents/A_nidulans/genome.size"
    
    # Generate tiles
    tiles = make_tiles_for_summit(summit_file, seq_range, chromatin_length_file)
    
    # Normalize tiles with ChIP-seq data
    norm_matrix = normalization_tiles_for_summits(chip_seq_file, tiles, seq_range, chromatin_length_file)
    
    # Plot the heatmap
    plot_heatmap(norm_matrix, max_density, 301)

# Example usage
path_folder = '/mnt/c/Users/59798/Documents/BAM_BED_files/ChIPseqused/'
path_folder_summit = '/mnt/c/Users/59798/Documents/BAM_BED_files/ChIPseqused/used_summits/'

acuKmyc_Ace_1 = f'{path_folder}AcuKmyc_Pro_ChIP_seq_rep1.bam.bed'
acuMHA_Ace_1 = f'{path_folder}AcuMHA_Pro_ChIP_seq_rep1.bam.bed'

KM_overlap = f'{path_folder_summit}AcuK_AcuM_summits_overlap.bed'
K_specific = f'{path_folder_summit}AcuK_AcuM_summits_AcuK_Specific.bed'
M_specific = f'{path_folder_summit}AcuK_AcuM_summits_AcuM_Specific.bed'

# Generate plots
plot_figures(KM_overlap, acuMHA_Ace_1, 600)
plot_figures(KM_overlap, acuKmyc_Ace_1, 600)
plot_figures(K_specific, acuMHA_Ace_1, 600)
plot_figures(K_specific, acuKmyc_Ace_1, 600)
plot_figures(M_specific, acuMHA_Ace_1, 600)
plot_figures(M_specific, acuKmyc_Ace_1, 600)
