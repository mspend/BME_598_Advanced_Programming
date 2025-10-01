# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pandas as pd

panCancer_phenos = pd.read_csv('phenotypes_panCancer.csv', header = 0, index_col=0)

print(panCancer_phenos.shape)

# Select all STAD tumors
panCancer_phenos_stad = panCancer_phenos.loc[panCancer_phenos['tumor'] == 'STAD',: ]

print(panCancer_phenos_stad.shape)


# Explor which columns have the most null values
print(panCancer_phenos_stad.info())
print(panCancer_phenos_stad.isnull().sum())

# Drop the columns we don't need to keep the most data possible
panCancer_phenos_stad = panCancer_phenos_stad.drop(columns = 'subtype')
print(panCancer_phenos_stad.isnull().sum())

print(panCancer_phenos_stad.shape)

# Drop missing values
panCancer_phenos_stad = panCancer_phenos_stad.dropna()

print(panCancer_phenos_stad.shape)



# Replace missing values with zeros
panCancer_phenos_noNA = panCancer_phenos.fillna(0)


# Load up leukocyte data
leukocyte_df = pd.read_csv('leukocyte_TCGA.csv', header = 0, index_col=0)


# Notice this DataFrame has the same number of rows as panCancer_phenos
print(leukocyte_df.shape)


# Merging two datasets
panCancer_phenos_merged = pd.merge(panCancer_phenos,leukocyte_df, on = 'bcr_patient_barcode')

# How many non-null values for TotalLeukocyte?
print(panCancer_phenos_merged.info())


# A groupby operation involves some combination of splitting the object, applying a function, and combining the results
# Note that groupby doesn't return anything unless it's applied to a function
median_leuk_by_tumor = panCancer_phenos_merged.groupby('tumor')['TotalLeukocyte'].median()
print(median_leuk_by_tumor.sort_values())




# Drop the columns we don't need to keep the most data possible
print(panCancer_phenos_merged.isnull().sum())
panCancer_phenos_filtered = panCancer_phenos_merged.drop(columns = 'subtype')
panCancer_phenos_filtered = panCancer_phenos_filtered.drop(columns = 'gender')
panCancer_phenos_filtered = panCancer_phenos_filtered.drop(columns = 'race')
panCancer_phenos_filtered = panCancer_phenos_filtered.drop(columns = 'OS')
panCancer_phenos_filtered = panCancer_phenos_filtered.drop(columns = 'OS.time')


print(panCancer_phenos_filtered.shape)

# Drop missing values
panCancer_phenos_filtered = panCancer_phenos_filtered.dropna()

print(panCancer_phenos_filtered.shape)
print(panCancer_phenos_filtered.info())






                     
                     
                     
                     
    