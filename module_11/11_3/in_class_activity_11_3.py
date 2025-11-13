## Import libraries
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
import scipy.stats as stats
import sklearn 
import GEOparse 

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.cm as cm
import matplotlib.pyplot as plt

# these lines are important so the text on the PDF doesn't get cut up in weird ways
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42

# clustering packages
from sklearn.preprocessing import StandardScaler, QuantileTransformer
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, silhouette_samples


# not sure if we need these
# from matplotlib.patches import Patch

# load up the data
gse = GEOparse.get_GEO('GSE11292')


## Get a list of GSMs
print(gse.gsms.keys())
print(f'The number of GSMs in this GSE is {len((gse.gsms.keys()))}.')

## Create an expression matrix
# Note: This won't work well with lots of samples!
expr = gse.pivot_samples('VALUE')
print(expr)

## Print out phenotype data
# Be sure to add the code below for colorizing based on time points later
print(gse.phenotype_data)
print(gse.phenotype_data.loc['GSM285027'])
convert_GSMs = gse.phenotype_data['title'].to_dict()
convert_GSMs = {i:convert_GSMs[i].split('_')[-1] for i in convert_GSMs} 

# Getting rid of the ThGFP and ThGARP samples
drop_me = ['GSM2850'+str(i) for i in range(45,51)]
expr2 = expr.drop(drop_me, axis=1)
print(f'After dropping the GFP and GARP, the dimensions of our expression matrix are {expr2.shape}.')
num_dropped = expr.shape[1] - expr2.shape[1]
print(f'Meaning we dropped {num_dropped} GSM samples.')

# log2 transform our data
# The data is very spread out
# if you don't log tansform the data, the data will be hard to plot
expr3 = np.log2(expr2)

# Make a boxplot of the log2 transformed data
with PdfPages('boxplot_GSE11292_pre_transform.pdf') as pdf:
    plt.boxplot(expr3)
    plt.xlabel('Samples')
    plt.ylabel('Expression (Log2(Signal))')
    pdf.savefig()

# Quantile normalize the data
expr4 = QuantileTransformer().fit_transform(expr3)
# Unfortunately, the above function converts expr3 into a numpy array, so we lose the index and columns.
# Convert it back to a DataFrame and add the index and columns
expr4 = pd.DataFrame(expr4)
expr4.index = expr3.index
expr4.columns = expr3.columns

# Make a boxplot of the quantile normalized data
with PdfPages('boxplot_GSE11292_post_transform.pdf') as pdf:
    plt.boxplot(expr4)
    plt.xlabel('Samples')
    plt.ylabel('Expression (Log2(Signal))')
    pdf.savefig()

## Feature selection
top3000 = expr4.var(axis=1).sort_values(ascending=False).index[range(3000)]

## Standardize data
tmp = StandardScaler().fit_transform(expr4.loc[top3000])

## Cluster using kMeans
# Compute silhouette scores and plots for k-means clustering applied across a range of 2 to 20 clusters
sil_km = []
with PdfPages('km_silhouettes_GSE11292.pdf') as pdf:
    for i in range(2,21):
        n_clusters = i
        km1 = KMeans(n_clusters=i).fit(tmp)

        # The silhouette_score gives the average value for all the samples.
        # This gives a perspective into the density and separation of the formed
        # clusters
        sil_km.append(silhouette_score(tmp, km1.labels_))

        # Create a subplot with 1 row and 2 columns
        fig, ax1 = plt.subplots(1, 1)
        #fig.set_size_inches(7, 7)

        # The silhouette coefficient can range from -1, 1
        ax1.set_xlim([-1, 1])

        # The (n_clusters+1)*10 is for inserting blank space between silhouette
        # plots of individual clusters, to demarcate them clearly.
        ax1.set_ylim([0, len(tmp) + (n_clusters + 1) * 10])

        # Compute the silhouette scores for each sample
        sample_silhouette_values = silhouette_samples(tmp, km1.labels_)

        y_lower = 10
        for j in range(i):
            # Aggregate the silhouette scores for samples belonging to
            # cluster i, and sort them
            jth_cluster_silhouette_values = \
                sample_silhouette_values[km1.labels_ == j]

            jth_cluster_silhouette_values.sort()

            size_cluster_j = jth_cluster_silhouette_values.shape[0]
            y_upper = y_lower + size_cluster_j

            color = cm.nipy_spectral(float(j) / n_clusters)
            ax1.fill_betweenx(np.arange(y_lower, y_upper),
                              0, jth_cluster_silhouette_values,
                              facecolor=color, edgecolor=color, alpha=0.7)

            # Label the silhouette plots with their cluster numbers at the middle
            ax1.text(-0.05, y_lower + 0.5 * size_cluster_j, str(j))

            # Compute the new y_lower for next plot
            y_lower = y_upper + 10  # 10 for the 0 samples

        ax1.set_title("The silhouette plot for the various clusters.")
        ax1.set_xlabel("The silhouette coefficient values")
        ax1.set_ylabel("Cluster label")

        # The vertical line for average silhouette score of all the values
        ax1.axvline(x=sil_km[i-2], color="red", linestyle="--")

        ax1.set_yticks([])  # Clear the yaxis labels / ticks
        ax1.set_xticks([-1,-0.8,-0.6,-0.4,-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1])

        # Save figure to pdf
        pdf.savefig(fig)
        plt.close()

    # Finally, plot silhouette score vs. number of clusters
    fig = plt.figure()
    plt.plot(range(2,21),sil_km)
    plt.xticks(range(2,21))
    plt.xlabel('Number of clusters')
    plt.ylabel('Average silhouette score')
    plt.title('Silhouette analysis for KMeans')
    pdf.savefig(fig)
    plt.close()


## Return an eigengene for a gene expression data given a set of genes
def getEigengene(gexp, genes):
    pca = PCA(n_components=1)
    gexp_pca = pd.DataFrame(pca.fit(gexp.loc[genes].T).transform(gexp.loc[genes].T), index = gexp.columns)
    eigengene = gexp_pca[0]
    if sum([stats.pearsonr(gexp.loc[i],eigengene)[0] for i in genes])/len(genes) > 0:
        return eigengene
    else:
        return -eigengene

# Run KMeans clustering
k = 6
km1 = KMeans(n_clusters=k, random_state=42).fit(tmp)

# Compute eigengenes for each cluster
eigengenes = pd.concat(
    [getEigengene(expr2.loc[top3000], top3000[km1.labels_ == i]) for i in range(k)],
    axis=1
)
eigengenes.columns = [f"Cluster {i+1}" for i in range(k)]

# Z-score eigengenes across samples
eigengenes_z = pd.DataFrame(
    StandardScaler().fit_transform(eigengenes),
    index=eigengenes.index,
    columns=eigengenes.columns
)

# Clip to highlight relative differences
eigengenes_z_clipped = eigengenes_z.clip(-2, 2)

# Define timepoint colors
time_labels = ['T' + str(i) + 'min' for i in range(0, 380, 20)]
colors = dict(zip(time_labels, sns.color_palette('coolwarm', n_colors=len(time_labels))))
col_colors = [colors[convert_GSMs[i]] for i in expr2.columns]

# Create clustermap
g = sns.clustermap(
    eigengenes_z_clipped.T,
    cmap="RdBu_r",
    col_colors=col_colors,
    col_cluster=False,
    center=0,
    vmin=-2, vmax=2,
    figsize=(10, 8)  # taller figure to fit legend + title comfortably
)

# Add a title BELOW the legend (using ax_heatmap title instead of suptitle)
g.ax_heatmap.set_title("Z-scored Eigengenes of K-Means Clusters (k=6)", fontsize=14, pad=30)

# Label axes
g.ax_heatmap.set_xlabel("Samples / Time Points", fontsize=12)
g.ax_heatmap.set_ylabel("Cluster Eigengenes", fontsize=12)

# Add legend for timepoints
for label in time_labels:
    g.ax_col_dendrogram.bar(0, 0, color=colors[label], label=label, linewidth=0)

# Position legend slightly higher so it doesn't overlap the title
g.ax_col_dendrogram.legend(
    loc="upper center",
    ncol=4,
    bbox_to_anchor=(0.5, 1.3),  # move legend higher
    title="Time Points",
    fontsize=8,
    title_fontsize=9
)

# Adjust layout so nothing overlaps
plt.subplots_adjust(top=0.9, bottom=0.1)

# Save to PDF
with PdfPages('eigengenes_GSE11292_highcontrast_tidy.pdf') as pdf:
    pdf.savefig(g.fig, bbox_inches='tight')
    plt.close(g.fig)