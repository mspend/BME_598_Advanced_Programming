#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 09:33:10 2025

@author: maurispendlove
"""

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

fev = pd.read_table(filepath_or_buffer='fev_dat.txt', delimiter='\t')
print(fev.size)
print(fev.shape)

fev_continuous_vars = fev.drop(columns=['sex','smoke'])

with PdfPages('fev_pairplots.pdf') as savedPDF:
    # Pairplot
    sns.pairplot(fev_continuous_vars)
    savedPDF.savefig()
    plt.close()


#### Linear regression ####

### Model 1: Does smoking affect FEV? ###

# Build the model
model1 = smf.ols('FEV ~ smoke', data = fev)

# Fit the model
results1 = model1.fit()

# Print results
print(results1.summary())

### Model 2: Effect of smoking and age on FEV ###

# Build the model
model2 = smf.ols('FEV ~ age + smoke', data = fev)

# Fit the model
results2 = model2.fit()

# Print results
print(results2.summary())


### Model 3: Effect of smoking, age, and height on FEV with height modeled as quadratic function ###

# Build the model
model3 = smf.ols('FEV ~ age + ht + np.power(ht, 2) + smoke', data = fev)

# Fit the model
results3 = model3.fit()

# Print results
print(results3.summary())


### Model 4: Effect of smoking, age, and their interaction on FEV ###

# Build the model
# Use Patsy's I() to include the squared term directly. This specifies smoke and age^2.
model4 = smf.ols('FEV ~ age*smoke', data = fev)

# Fit the model
results4 = model4.fit()

# Print results
print(results4.summary())


### Model 5: Interaction between age and smoking, and then add quadratic height on FEV ###

# Build the model
model5 = smf.ols('FEV ~ age * smoke + ht + np.power(ht, 2)', data=fev)

# Fit the model
results5 = model5.fit()

# Print results
print(results5.summary())


# formulas I have attempted for model 3 that are listed as incorrect by autograder:

# model3 = smf.ols('FEV ~ smoke + I(age ** 2) + ht', data = fev)
# model3 = smf.ols('FEV ~ smoke + age + I(ht ** 2)', data = fev)
# model3 = smf.ols('FEV ~ smoke + age + ht + I(ht ** 2)', data = fev)
# model3 = smf.ols('FEV ~ smoke + age + ht', data = fev)
# model3 = smf.ols('FEV ~ smoke + age + ht*ht', data = fev)
# model3 = smf.ols('FEV ~ smoke + age + I(ht^2)', data = fev)
# model3 = smf.ols('FEV ~ smoke + age + ht:ht', data = fev)