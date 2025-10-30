#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 09:33:10 2025

@author: maurispendlove
"""

import matplotlib
# matplotlib.use('Agg')
# matplotlib.rcParams['pdf.fonttype'] = 42
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

fev = pd.read_table(filepath_or_buffer='fev_dat.txt', delimiter='\t')
print(fev.size)

fev_continuous_vars = fev.drop(columns=['sex','smoke'])

with PdfPages('fev_pairplots.pdf') as savedPDF:
    # Pairplot
    sns.pairplot(fev_continuous_vars)
    savedPDF.savefig()
    plt.close()


## Linear regression

# # Build the model
# model1 = smf.ols('FEV ~ ', data = fev)

# # Fit the model
# results1 = model1.fit()

# # Collect results
# print(results1.summary())









# look for an interaction between age and smoking
# Question 7
# smf.ols('FEV~age*smoke',data=fev)