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

# What disease states are present in the training dataset?
print(phenosTrain['disease_state'].value_counts())

# Get rid of lung cancer and Baseline 
# (How is baseline different from control? We don't know)
subset_train = phenosTrain.index[phenosTrain['disease_state'].isin(['Control','Active Sarcoid','TB','Non-active sarcoidosis'])]
phenosTrain = phenosTrain.loc[subset_train]

# Number of each disease_state in GSE42830
print(phenosTrain['disease_state'].value_counts())

# Phenotypes for test dataset (GSE42826)
phenosTest = gseTest.phenotype_data[['characteristics_ch1.0.gender','characteristics_ch1.1.ethnicity','characteristics_ch1.2.disease state']]
phenosTest.columns = ['gender','ethnicity','disease_state']

# What disease states are present in the testing dataset?
print(phenosTest['disease_state'].value_counts())

# Get rid of cancer and pneumonia
subset_test = phenosTest.index[phenosTest['disease_state'].isin(['Control','Active Sarcoid','TB','Non-active sarcoidosis'])]
phenosTest = phenosTest.loc[subset_test]
print(phenosTest['disease_state'].value_counts())

# Phenotypes for validation dataset (GSE42825)
phenosValidation = gseValidation.phenotype_data[['characteristics_ch1.0.gender','characteristics_ch1.1.ethnicity','characteristics_ch1.2.disease state']]
phenosValidation.columns = ['gender','ethnicity','disease_state']

# What disease states are present in the validation dataset?
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

## 5. Feature Selection
top1000 = gexpTrain.var(axis=1).sort_values(ascending=False).index[range(1000)]

## 6. Preprocess: standard scaling (important for NN on tabular features)
# Why do you use fit_transform on the train data and transform on the test data?
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(gexpTrain.loc[top1000].T)
X_test_scaled  = scaler.transform(gexpTest.loc[top1000].T)

# Define the y values
y_tran = 



# Define the y values
convertMe = {'TB': 0, 'Active Sarcoid': 1, 'Non-actives sarcoidosis': 2, 'Control': 3}
y_train = to_categorical([convertMe[i] for i in phenosTrain['disease state']],4])
y_test = to_categorical([convertMe[i] for i in phenosTest['disease state']],4])

## 7. Training the model
history = model.fit(
    X_train_scaled, y_train,
    validation_split=0.15,
    epochs=100,
    batch_size=32,
    verbose=2
)

## 8. Evaluate on test set: metrics important in biomedical context
y_proba = model.predict(X_test_scaled).ravel()
y_pred = (y_proba >= 0.5).astype(int)

acc = accuracy_score(y_test, y_pred)
prec = precision_score(y_test, y_pred)
recall = recall_score(y_test, y_pred)   # sensitivity / recall
auc = roc_auc_score(y_test, y_proba)

print(f"Test accuracy: {acc:.4f}")
print(f"Test precision: {prec:.4f}")
print(f"Test recall (sensitivity): {recall:.4f}")
print(f"Test ROC AUC: {auc:.4f}")

# Confusion matrix: [ [TN, FP], [FN, TP] ]
cm = confusion_matrix(y_test, y_pred)
tn, fp, fn, tp = cm.ravel()
specificity = tn / (tn + fp)
print("Confusion matrix:\n", cm)
print(f"Specificity: {specificity:.4f}")

## 9. Plot ROC and training curves
with PdfPages('training_curves_multiclass.pdf') as pdf:
    fig, ax = plt.subplots(1,2, figsize=(12,4))
    ax[0].plot(history.history['loss'], label='train_loss')
    ax[0].plot(history.history['val_loss'], label='val_loss')
    ax[0].set_xlabel('Epoch')
    ax[0].set_ylabel('Loss')
    ax[0].legend()
        
    if 'accuracy' in history.history:
        ax[1].plot(history.history['accuracy'], label='train_acc')
        ax[1].plot(history.history['val_accuracy'], label='val_acc')
        ax[1].set_xlabel('Epoch')
        ax[1].set_ylabel('Accuracy')    
        ax[1].legend()
    pdf.savefig()
    plt.close()
    
with PdfPages('ROC_curve_multiclass.pdf') as pdf:
    # ROC
    fpr, tpr, _ = roc_curve(y_test, y_proba)
    plt.figure(figsize=(6,6))
    plt.plot(fpr, tpr, label=f"AUC = {roc_auc_score(y_test, y_proba):.3f}")
    plt.plot([0,1], [0,1], linestyle='--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate (Sensitivity)')
    plt.title('ROC Curve')
    plt.legend()
    pdf.savefig()
    plt.close()


## 10. Save model in keras format
model.save('breast_cancer_nn_savedmodel.keras')

