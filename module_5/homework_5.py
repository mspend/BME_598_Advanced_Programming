import pandas as pd
import matplotlib.pyplot as plt

# What is the  impact of age on a set of genes in 
# bladder urothelial carcinoma, lung adenocarcinoma, and skin cutaneous melanoma?

# Are there sex-specific effects (gender) on overall survival (OS.time)

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

