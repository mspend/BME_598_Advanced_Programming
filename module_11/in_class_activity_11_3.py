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
from sklearn import cluster, datasets, mixture
from sklearn.neighbors import kneighbors_graph
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler, QuantileTransformer
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, SpectralClustering, AgglomerativeClustering
from sklearn.mixture import GaussianMixture
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
    pdf.close()

# Quantile normalize the data
expr4 = QuantileTransformer().fit_transform(expr3)
# Unfortunatley, the above function converts expr3 into a numpy array, so we lose the index and columns.
# Convert it back to a DataFrame and add the index and columns
expr4 = pd.DataFrame(expr4)
expr4.index = expr3.index
expr4.columns = expr3.columns

# Make a boxplot of the log2 transformed data
with PdfPages('boxplot_GSE11292_post_transform.pdf') as pdf:
    plt.boxplot(expr4)
    plt.xlabel('Samples')
    plt.ylabel('Expression (Log2(Signal))')
    pdf.savefig()
    pdf.close()


## Feature selection
top3000 = expr4.var(axis=1).sort_values(ascending=False).index[range(3000)]












## Scaling
tmp = StandardScaler().fit_transform(expr2.loc[top3000])


# dictionry comprehension
convert_GSMs = gse.phenotype_data['title'].to_dict()
# grab the last element of the title
convert_GSMs = {i:convert_GSMs[i].split('_')[-1] for i in convert_GSMs}
