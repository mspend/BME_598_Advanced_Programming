#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 09:39:06 2025

@author: maurispendlove
"""


# Load up libraries
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

# Set up the plot
fig, ax = plt.subplots(1,1)
# plt.show()

# Creates a curve based on PDF
x = np.linspace(norm.ppf(0.01), norm.ppf(0.99), 100)
ax.plot(x, norm.pdf(x), 'r-', lw=5, alpha=0.6, label='empirical')
# plt.show()

# Alternative method by directly creating a normal distribution instance
rv = norm()
ax.plot(x, rv.pdf(x), 'k-', lw=2, label='direct')
# plt.show()

# Make random variable
norm_thousand = norm.rvs(size=1000)
ax.hist(norm_thousand, density=True, bins='auto', histtype='stepfilled',alpha=0.2)
ax.set_xlim([x[0], x[-1]])
ax.legend(loc='best', frameon=False)
plt.show()
# plt.savefig('normal_distribution.pdf')




# generates a single random variable / sample
stats.norm.rvs()

# generates an array of 5 random variables
stats.norm.rvs(size=5)

stats.norm.pdf(x, loc=0, scale=1)


# Multipage PDF

from matplotlib.backends.backend_pdf import PdfPages

# Create a PdfPages object
with PdfPages('my_multipage_plots.pdf') as pdf:
    
    # First plot
    plt.figure()
    plt.plot([1, 2, 3], [4, 5, 6])
    plt.title('Plot One')
    pdf.attach_note('Plot One')
    pdf.savefig() # Save the current figure to a new page
    plt.close() # Close the figure to free up memory
    
    # Second plot
    plt.figure()
    plt.bar(['A', 'B', 'C'], [10, 20, 15])
    plt.title('Plot Two')
    pdf.savefig() # Save the current figure to a new page
    plt.close() # Close the figure




