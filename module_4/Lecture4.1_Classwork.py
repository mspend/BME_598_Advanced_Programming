#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  9 09:27:09 2025

@author: maurispendlove
"""

import pandas as pd

food_df = pd.DataFrame({'Food': {'Bill': 'pizza', 'Joe':
'beans', 'Umu': 'hamburger'},'Drink':{'Bill':'jolt_cola','Joe': 'tea', 'Umu': 'coca_cola' }})
    
food_df.to_csv('/Users/maurispendlove/ASU Dropbox/Mauri Spendlove/Classes/BME 598 Programming/Data/lunch_DataFrame.csv')


panCancer_phenos = pd.read_csv('/Users/maurispendlove/ASU Dropbox/Mauri Spendlove/Classes/BME 598 Programming/Data/phenotypes_panCancer.csv',index_col='bcr_patient_barcode')


print(panCancer_phenos.info())

# descriptive statistics
print(panCancer_phenos.describe())

# How many unique tumor types are there
print(panCancer_phenos['tumor'].unique())

# How many counts are there of each tumor type, descending?

print(panCancer_phenos['tumor'].value_counts(sort=True, ascending=False, dropna=True))
# value_counts returns a Series containing the frequency of each distinct row in the Dataframe
# the 3 parameters are unecessary because that's their default value

# This also works
print(panCancer_phenos.value_counts(subset='tumor'))


# Percent of GBM patients that are male
percent = len(panCancer_phenos.loc[panCancer_phenos.tumor=='GBM',:][panCancer_phenos.gender=='MALE']) /len(panCancer_phenos.loc[panCancer_phenos.tumor == 'GBM',:])
print(percent)







# First five rows
print(panCancer_phenos.head())
# Last five rows
print(panCancer_phenos.tail())

# print(panCancer_phenos.loc[:, 'tumor'].value_counts())
                                             
# panCancer_phenos_stad_yt40 = panCancer_phenos.loc[(panCancer_phenos['tumor'] == 'STAD’) & (panCancer_phenos['age_at_initial_pathologic_diagnosis'] <= 40), :]                                            
                                                   
# panCancer_phenos.loc[panCancer_phenos.loc[:,'tumor'].isin(['ESCA','STAD','COAD','READ']),:]

# panCancer_phenos.loc[panCancer_phenos.loc[:,'subtype'].isin(['MSI', 'CIN'])]

