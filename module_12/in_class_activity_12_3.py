##########################################################
## OncoMerge:  in_class_activity_12_3.py                ##
##  ______     ______     __  __                        ##
## /\  __ \   /\  ___\   /\ \/\ \                       ##
## \ \  __ \  \ \___  \  \ \ \_\ \                      ##
##  \ \_\ \_\  \/\_____\  \ \_____\                     ##
##   \/_/\/_/   \/_____/   \/_____/                     ##
## @Developed by: Plaisier Lab                          ##
##   (https://plaisierlab.engineering.asu.edu/)         ##
##   Arizona State University                           ##
##   242 ISTB1, 550 E Orange St                         ##
##   Tempe, AZ  85281                                   ##
## @Author:  Chris Plaisier                             ##
## @License:  GNU GPLv3                                 ##
##                                                      ##
## If this program is used in your analysis please      ##
## mention who built it. Thanks. :-)                    ##
##########################################################


## May need to install
# pip install GEOparse
# pip install tensorflow


## 1. Load libraries
import GEOparse as gp
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

import tensorflow as tf
from tensorflow.keras.utils import to_categorical
from tensorflow.keras import layers, models, optimizers
from sklearn.metrics import (
    roc_auc_score, roc_curve, confusion_matrix,
    accuracy_score, precision_score, recall_score
)


####################################################################
## Load up real-world data: Human whole blood microarray study    ##
## to compare patients with tuberculosis, sarcoidosis,            ##
## pneumonia, and lung cancer                                     ##
## PMID = 23940611                                                ##
## Note:  The authors have split their data into train, test,     ##
## and validation. But as we will learn using cross-validation    ##
## is a better approach. We will try both methods.                ##
####################################################################

## 2. Load up real-world data into pandas
# Using data from super-series GSE42834 on GEO:
#  - Training = https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42830
#  - Test = https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42826
#  - Validation = https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42825
gseTrain = gp.get_GEO(filepath="data/GSE42830_family.soft.gz")
gseTest = gp.get_GEO(filepath="data/GSE42826_family.soft.gz")
gseValidation = gp.get_GEO(filepath="data/GSE42825_family.soft.gz")


## 3. Extract phenotypes/labels (disease state in this case are the labels)
## For starters let's include four disease states:  Control, Active Sarcoid, TB, and Non-active sarcoidosis
print(gseTrain.phenotype_data)

# Select out gender, ethnicity, and disease state
phenosTrain = gseTrain.phenotype_data[['characteristics_ch1.0.gender','characteristics_ch1.1.ethnicity','characteristics_ch1.2.disease state']]
phenosTrain.columns = ['gender','ethnicity','disease_state']
phenosTrain

# Get rid of cacner and 
subset_train = phenosTrain.index[phenosTrain['disease_state'].isin(['Control','Active Sarcoid','TB','Non-active sarcoidosis'])]
phenosTrain = phenosTrain.loc[subset_train]

# Number of each disease_state in GSE42830
print(phenosTrain['disease_state'].value_counts())

# Phenotypes for test dataset (GSE42826)
phenosTest = gseTest.phenotype_data[['characteristics_ch1.0.gender','characteristics_ch1.1.ethnicity','characteristics_ch1.2.disease state']]
phenosTest.columns = ['gender','ethnicity','disease_state']
print(phenosTest['disease_state'].value_counts())

# Get rid of cacner and 
subset_test = phenosTest.index[phenosTest['disease_state'].isin(['Control','Active Sarcoid','TB','Non-active sarcoidosis'])]
phenosTest = phenosTest.loc[subset_test]

# Phenotypes for validation dataset (GSE42825)
phenosValidation = gseValidation.phenotype_data[['characteristics_ch1.0.gender','characteristics_ch1.1.ethnicity','characteristics_ch1.2.disease state']]
phenosValidation.columns = ['gender','ethnicity','disease_state']
print(phenosValidation['disease_state'].value_counts())
# Uh-oh! Active sarcoidosis does not equal Active Sarcoid, from GSE42830 and GSE42826
# Let's fix it by harmonizing the validation to Active Sarcoid
phenosValidation.loc[phenosValidation.disease_state=='Active sarcoidosis','disease_state'] = 'Active Sarcoid'
print(phenosValidation['disease_state'].value_counts())


## 4. Extract gene expression data
# What columns are available per sample
# print(gseTrain.gsms['GSM1050928'].columns)
# We want to have VALUE for each sample combined into one dataframe

# Build a DataFrame from VALUEs using the pivot_samples function
# With genes as the rows and samples as the columns
gexpTrain = gseTrain.pivot_samples('VALUE').loc[:,subset_train]
gexpTest = gseTest.pivot_samples('VALUE').loc[:,subset_test]
gexpValidation = gseValidation.pivot_samples('VALUE')

