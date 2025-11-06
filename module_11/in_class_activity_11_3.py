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

# these lines are important so the text on the PDF doesn't get cut up in weird ways
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42

# clustering packages
from sklearn import cluster, datasets, mixture
from sklearn.neighbors import kneighbors_graph
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, SpectralClustering, AgglomerativeClustering
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score, silhouette_samples
from scipy import stats



# not sure if we need these
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

