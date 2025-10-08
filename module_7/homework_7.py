#############################################################
## BME 598:  in_class_activity_7_2.py                      ##
##  ______     ______     __  __                           ##
## /\  __ \   /\  ___\   /\ \/\ \                          ##
## \ \  __ \  \ \___  \  \ \ \_\ \                         ##
##  \ \_\ \_\  \/\_____\  \ \_____\                        ##
##   \/_/\/_/   \/_____/   \/_____/                        ##
## @Developed by: Plaisier Lab                             ##
##   (https://plaisierlab.engineering.asu.edu/)            ##
##   Arizona State University                              ##
##   242 ISTB1, 550 E Orange St                            ##
##   Tempe, AZ  85281                                      ##
## @Author:  Chris Plaisier                                ##
## @License:  GNU GPLv3                                    ##
##                                                         ##
## If this program is used in your analysis please         ##
## mention who built it. Thanks. :-)                       ##
#############################################################

## Total points = 100 pts (70 pts through autograder, 30 pts through PDF submitted through Canvas)

## Your name
# Mauri Spendlove

## 1. (10pts) Prepare an executive summary that clearly states the goals of the proposed analysis and provides a concise overview of the analytical steps to be undertaken�
#             A. The executive summary needs to be written in plain text; not as bullet points. And should be no more than one paragraph.
#             B. The overview of analytical steps should be written as a outline with bullet points. Provide more detail than the headings of the Deliverables section.


## 2. (10pts) Load up the packages you will need to run
#  Deliverables:
#     a. Restart console: Consoles tab -> Restart Kernel
#     b. Import all the packages
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


## 3. (10pts) Load up and inspect data and function
#  Deliverables:
#     a. Use a text editor to determine the structure of the 'Raw_data_TRP_with_ischemic_stroke.csv' file
#        - Delimiter?
#        - Header?
#        - Index column?
#     c. Read the file using pandas into a variable named 'stroke'
#     d. Make sure the shape of the file is correct
#           Expectation is (519 rows, by 15 columns)
stroke = pd.read_csv('Raw_data_TRP_with_ischemic_stroke.csv', delimiter=',', header=0)
print(stroke.shape) 

## 4. (10pts) Make a pairs plot of your data to identify which variables to test in correlation studies.
#  Deliverables:
#     a. Write a pairsplot PDF using PdfPages to a pdf called 'stroke_data_pairplot.pdf'
with PdfPages('stroke_data_pairplot.pdf') as pdf:
    plt.figure(figsize=(6, 6))
    sns.pairplot(data=stroke)
    plt.tight_layout()
    pdf.savefig()
    plt.close()

## 5. (20pts) Compute pairwise correlations
#  Deliverables:
#     a. Make a variable 'quant_variables' that is a list of all the columns that should be inlcuded in
#        the corrleation studies
#     b. Make two dictionaries to hold the results of correlation analysis:
#        - 'pearson_cor' to hold the correlation values
#        - 'pearsons_pv' to hold the p-values
#     c. Fill the dictionaries up with the values from the pairwise correlation comparisons of the 'quant_variables'
#        - NOTE: set the comparison of a variable to itself as 0

# Notice that the data has NA values. This will give an error when we run correlations later.
print(stroke.info())
stroke_filtered = stroke.dropna()

vars_as_string = "Age,SBP,DBP,HbA1c,FBS,TC,LDL,HDL,TG,WBC,CRP,AIP,TyG"

quant_variables = vars_as_string.split(sep=',')
pearson_cor = {}
pearson_pv = {}

for q1 in quant_variables:
    pearson_cor[q1] = {}
    pearson_pv[q1] = {}
    for q2 in quant_variables:
        if not q1==q2:
            pearson_cor[q1][q2], pearson_pv[q1][q2] = pearsonr(stroke_filtered[q1], stroke_filtered[q2])
        else:
            pearson_cor[q1][q2], pearson_pv[q1][q2] = (0, np.nan)


## 6. (10pts) Convert result dictionaries into pd.DataFrames and write them out as CSV files
#  Deliverables:
#     a. Convert both result dictionaries into pd.DataFrames:
#        - 'pearson_cor' -> 'pearson_cor_df'
#        - 'pearson_pv' -> 'pearson_pv_df'
#     b. Concatenate the two pd.DataFrames and change the names of the columns to denote correlation coefficient 'R',
#        and p-value '.pv'
#     c. Then write out the merged pd.DataFrame as 'stroke_pearson_correlation.csv'

# Convert to pandas DataFrames
pearson_cor_df = pd.DataFrame(pearson_cor)
pearson_pv_df = pd.DataFrame(pearson_pv)

# Confirm the conversion to DataFrame worked
print(type(pearson_cor_df))

# Write out combined R and PV matrix
# axis = 1 means concatenate along the columns (axis=0 is for index)
combined_df = pd.concat([pearson_cor_df, pearson_pv_df], axis=1).sort_index(axis=1)

# add labels to column headers to make the combined_df easier to read
# Otherwise there would be 2 columns with the same name for each variable
# .pv for P-value if the column numbers (indexed from 1) are even, .R for Pearson R
# When the dfs were concatenated, the df listed first (pearson_cor_df) was on the left, before sorting.
combined_df.columns = [combined_df.columns[i]+'.pv' if ((i+1)%2)==0 else combined_df.columns[i]+'.R' for i in range(len(combined_df.columns))]
combined_df.to_csv('stroke_pearson_correlation.csv')

## 7. (10pts) Visiualize the correlation matrix as a heatmap
#  Deliverables:
#     a. Generate a heatmap as a PDF called 'stroke_pearson_cor_df_heatmap.pdf'
#        - NOTE: Save the sns.heatmap to a variable called 'hm1' to facilitate the autograding, i.e.:
#                hm1 = sns.heatmap(<fill with your parameters!>)
#     b. Use the 'vlag' cmap which goes from red to white to blue
#     c. Maximum and minimum correlation values should be set at 1 and -1, respectively
#     d. Must include units for the colorbar legend

# Make a plot of the relationships
with PdfPages('stroke_pearson_cor_df_heatmap.pdf') as pdf:
    plt.figure(figsize=(6, 6))
    hm1 = sns.heatmap(data=pearson_cor_df, cmap=sns.color_palette("vlag", as_cmap=True), vmin=-1, vmax=1, cbar_kws={'label': 'Pearson\'s R'})
    plt.tight_layout()
    pdf.savefig()
    plt.close()

## 8. (20pts) Biological interpretation that is to be completed in a separate Word document that will be converted into a PDF for submission
#  Deliverables:
#     a. Using a minimum of 3 and a maximum of 5 paragraphs of plain text (not bullet points), describe how the data in the figures, and
#        output from statistical tests in the results tables answer the following questions:
#        - How many relationships were observed to be significant?
#            - How many were positive?
#            - How many were negative?
#        - What is the meaning of the relationships?
#        - Any trends that could reflect tumor biology?
#        - Are the relationships known? Do some research to validate at least two relationships using one citation.
#     c. Please embed the figures and excerpts of the table to reinforce your interpretation.




