# Mauri Spendlove
# 

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

##### OBJECTIVES #####
# What is the  impact of age on a set of genes in 
# bladder urothelial carcinoma, lung adenocarcinoma, and skin cutaneous melanoma?

# Are there sex-specific effects (gender) on overall survival (OS.time)


# load up data
panCancer_phenos = pd.read_csv('phenotypes_panCancer.csv')

# first 5 rows
print(panCancer_phenos.head())

# names of columns
print(panCancer_phenos.columns)

# count of non-null data per column
print(panCancer_phenos.info())

# descriptive statistics
print(panCancer_phenos.describe())

# key variables for our analysis:  tumor, age_at_initial_pathologic_diagnosis, gender, and OS.time

# How many unique tumor types are there
print(panCancer_phenos['tumor'].unique())

# How many counts are there of each tumor type, descending?
print(panCancer_phenos.value_counts(subset='tumor'))

# Drop the columns we don't need to keep the most data possible
panCancer_phenos = panCancer_phenos.drop(columns = 'subtype')
panCancer_phenos = panCancer_phenos.drop(columns = 'OS')


# Cancers of interest: 
# 1. bladder urothelial carcinoma = BLCA
# 2. lung adenocarcinoma = LUAD
# 3. skin cutaneous melanoma = SKCM

cancers_of_interest = panCancer_phenos.loc[(panCancer_phenos.tumor=='BLCA') | (panCancer_phenos.tumor=='LUAD') | (panCancer_phenos.tumor=='SKCM'),:]

# Drop missing values
cancers_of_interest = cancers_of_interest.dropna()

##### Plot #1: Tumor Type vs Age
# Select only cancers of interest
cancers_of_interest = ["BLCA", "LUAD", "SKCM"]
subset = panCancer_phenos[panCancer_phenos["tumor"].isin(cancers_of_interest)]

# Prepare data for age violinplot
data_1 = [subset[subset["tumor"] == g]["age_at_initial_pathologic_diagnosis"].dropna()
        for g in cancers_of_interest]

# Create violin plot of age vs tumor type
plt.figure(figsize=(6,5))
plt.violinplot(data_1, showmeans=True)

# Label axes
plt.xticks(np.arange(1,len(cancers_of_interest)+1), cancers_of_interest)
plt.xlabel("Tumor Type")
plt.ylabel("Age at Initial Diagnosis (Years)")
plt.title("Examining Age at Initial Diagnosis")
plt.show()

##### Plot #2: Sex vs Survival
# Select data
sexes = ["MALE", "FEMALE"]
subset_2 = panCancer_phenos[panCancer_phenos["gender"].isin(sexes)]

# Prepare data for gender violinplot
data_2 = [panCancer_phenos[panCancer_phenos["gender"] == g]["OS.time"].dropna()
        for g in sexes]

# Create violin plot of gender vs survival time
plt.figure(figsize=(6,5))
plt.violinplot(data_2, showmeans=True)

# Label axes
plt.xticks(np.arange(1,len(sexes)+1), sexes)
plt.xlabel("Sex")
plt.ylabel("Survival (Days)")
plt.title("Relationship between Sex and Overall Survival time")
plt.show()

# Print median survival time for both men and women
median_survival = panCancer_phenos.groupby('gender')['OS.time'].median()
print(median_survival)

##### Plot #3: Race vs Age
# Select data
races = ['WHITE', 'BLACK OR AFRICAN AMERICAN', 'ASIAN', 
         'AMERICAN INDIAN OR ALASKA NATIVE', 'NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER', '[Unknown]']
subset_3 = panCancer_phenos[panCancer_phenos["race"].isin(races)]

# Prepare data: group age values by race
data_3 = [
    subset_3[subset_3["race"] == race]["age_at_initial_pathologic_diagnosis"].dropna()
    for race in races
]

# Create violin plot of age vs race
plt.figure(figsize=(6,5))
plt.violinplot(data_3, vert=False, showmeans=True)

# Label axes
plt.yticks(np.arange(1,len(races)+1), races)
plt.ylabel("Race")
plt.xlabel("Age at Initial Diagnosis (Years)")
plt.title("Relationship Between Race and Age at Initial Diagnosis")
plt.show()