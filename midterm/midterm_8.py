#############################################################
## BME 598:  midterm_8.py                                  ##
##  ______     ______     __  __                           ##
## /\  __ \   /\  ___\   /\ \/\ \                          ##
## \ \  __ \  \ \___  \  \ \ \_\ \                         ##
##  \ \_\ \_\  \/\_____\  \ \_____\                        ##
##   \/_/\/_/   \/_____/   \/_____/                        ##
## @Developed by: Plaisier Lab                             ##
##   (https://plaisierlab.engineering.asu.edu/)            ##
##   Arizona State University                              ##
##   242 ISTB1, 550 E Orange St                            ##
##   Tempe, AZ  85281                                      ##
## @Author:  Chris Plaisier                                ##
## @License:  GNU GPLv3                                    ##
##                                                         ##
## If this program is used in your analysis please         ##
## mention who built it. Thanks. :-)                       ##
#############################################################

## Import
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, ttest_ind
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sns


## Load up data
# Original Data from: https://archive.ics.uci.edu/ml/datasets/Heart+Disease Changes
# header but no index
# default is to add an index
clev = pd.read_table("proc_heart_cleve_3_withheader.tab", delimiter='\t',header=0) # <<-- Add any necessary code to load up data correctly
print('clev.shape = '+str(clev.shape))


## Compute ttest p-value for ST_dep_by_exerc vs. Disease
# Disease = -1 means they don't have the disease (control) = b
# Disease = 1 means they have the disease (case) = a
qualitative_trait = 'Disease' # <<-- Replace with qualitative trait name
quantitative_trait = 'ST_dep_by_exerc' # <<-- Repalce with quantitative trait name
print(ttest_ind(clev.loc[clev[qualitative_trait]==1,quantitative_trait], clev.loc[clev[qualitative_trait]==-1,quantitative_trait]))


## Compute Pearson correalation between maximum heart rate and ST_dep_by_exerc
qualitative_trait_1 = 'Max_heart_rate' # <<-- Replace with qualitative trait name
quantitative_trait_2 = 'ST_dep_by_exerc' # <<-- Repalce with quantitative trait name
print(pearsonr(clev[qualitative_trait_1], clev[quantitative_trait_2]))