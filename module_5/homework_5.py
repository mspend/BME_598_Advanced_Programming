import pandas as pd
import matplotlib.pyplot as plt

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

# Cancers of interest: 
# 1. bladder urothelial carcinoma = BLCA
# 2. lung adenocarcinoma = LUAD
# 3. skin cutaneous melanoma = SKCM

cancers_of_interest = panCancer_phenos.loc[(panCancer_phenos.tumor=='BLCA') | (panCancer_phenos.tumor=='LUAD') | (panCancer_phenos.tumor=='SKCM'),:]

# Drop the columns we don't need to keep the most data possible
cancers_of_interest = cancers_of_interest.drop(columns = 'subtype')
cancers_of_interest = cancers_of_interest.drop(columns = 'OS')

# Drop missing values
cancers_of_interest = cancers_of_interest.dropna()

