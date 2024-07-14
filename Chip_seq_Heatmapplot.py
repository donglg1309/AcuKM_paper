import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def plot_heatmap_B(x, y, z, k, filename):
    M_0h_p_matrix2 = np.copy(x)
    M_0h_p_matrix2[M_0h_p_matrix2 > z] = z
    y_input = np.copy(y)
    y_input[y_input > z] = z
    
    # Order rows by the mean of y_input in ascending order
    ordered_indices = np.argsort(np.mean(y_input, axis=1))
    ordered_matrix = M_0h_p_matrix2[ordered_indices, :]

    plt.figure(figsize=(10, 8))
    ax = sns.heatmap(
        ordered_matrix.T,  # Transpose to match the image plot orientation
        cmap="YlOrBr",     # Yellow to orange color palette
        cbar_kws={'label': 'Density'},
        vmin=0,
        vmax=z
    )
    ax.invert_yaxis()  # Invert the y-axis to match the R plot orientation
    plt.xlabel('')
    plt.ylabel('')
    plt.xticks([])
    plt.yticks([])
    plt.savefig(filename)
    plt.close()

# Example usage
x = np.random.rand(100, 301)  # Example data, replace with actual data
y = np.random.rand(100, 301)  # Example data, replace with actual data
z = 600
k = 301
filename = 'name.png'

plot_heatmap_B(x, y, z, k, filename)
