##########################################################
## BME 494/598:  Homework 6 Python code                 ##
##  ______     ______     __  __                        ##
## /\  __ \   /\  ___\   /\ \/\ \                       ##
## \ \  __ \  \ \___  \  \ \ \_\ \                      ##
##  \ \_\ \_\  \/\_____\  \ \_____\                     ##
##   \/_/\/_/   \/_____/   \/_____/                     ##
## @Developed for: BME 494/598 Applied Programming      ##
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

## Total points = 100 pts



## Your name
# Mauri Spendlove


## 1. (20pts) Load up the packages you will need to run
#  Deliverables:
#     a. Restart console: Consoles tab -> Restart Kernel
#     b. Import all the packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats


## 2. (10pts) Load up and inspect data and function
#  Deliverables:
#     a. Use a text editor to determine the structure of the 'proliferation_v3_w_redos_clean.csv' file
#        - Delimiter?
#        - Header?
#        - Index column?
#     c. Read the file using pandas into a variable named 'prolif'
#     d. Make sure the shape of the file is correct
#        - Expectation is (42 rows, by 14 columns)

prolif = pd.read_csv('proliferation_v3_w_redos_clean.csv', delimiter=',', header=0, index_col=0)
print(prolif.shape) 

## 3. (15pts) Background subtraction:
#  Deliverables:
#         The results for these proliferation assays were conducted on different plates, and therefore
#     the No_cell values will vary across plates. Thus, we must compute these for the cell line and
#     plate on which it was run. The different names in the No_cell column of prolif correspond to the
#     different plates, and the specific cell line it was compute for are in the columns labeled by
#     the cell line replicates (*.1, *.2, *.3, *.4).
#     a. Compute the median 'No cell' background for each cell line and store in a nested dictionary
#        named 'median_no_cell' that has name of the 'No_cell' row name and then the name of each cell
#        line under that:
#            Example:
#
#                median_no_cell = {'No_cell': {'HA':1234, 'T98G': 1234, 'U251': 1234},
#                                   ...
#                                 }
#
#     b. Then, subtract the median 'No cell' background value from each perturbation from the same cell
#        line. Make a new pandas DataFrame with corrected cell counts called 'prolif_corrected'

# The background is represented by the No_cell rows in the DataFrame
median_no_cells = {}
celllines = ['HA', 'T98G', 'U251']

# Pull out the 3 rows that correspond to No_cell plate
for index in prolif.index:
    if index.startswith("No_Cell"):
        # Pull out the row associated with that no cell plate
        plate = prolif.loc[index]
        temp_dict = {}
        # Iterate through the 3 cell lines and the 4 replicates in each cell line
        for cellline in celllines:
            replicates = [cellline+'.'+str(i) for i in range(1,5)]
            # Make a list with the fluorescent cell numbers for each cell line for each no cell plate
            cell_numbers = []
            for replicate in replicates:
                cell_numbers.append(plate[replicate])
                # print(plate[replicate])
                median = np.median(cell_numbers)
                temp_dict[cellline] = median
            # print(temp_dict)
        median_no_cells[index] = temp_dict
            # print(cell_numbers)
print(median_no_cells) 

# Make a copy of the prolif dataframe on which we will perform background subtraction.
prolif_corrected = prolif.copy()

# Iterate through every line in the dataframe
for index in prolif_corrected.index:
    # We want to calculate background subtraction on everything except for the no cell controls
    if not index.startswith("No_Cell"):
        # Turns each line into a series
        row = prolif_corrected.loc[index]

        for cellline in celllines:
            for row_index in row.index:
                if row_index.startswith(cellline):
                    value = row[row_index]
                    corrected_value = value - median_no_cells[row['No_cell']][cellline]
                    # replace values the dataframe with corrected values
                    prolif_corrected[row_index][index] = corrected_value


## 4. (15pts) Plot proliferation per treatment and cell line
#  Deliverables:
#     a. Convert 'prolif_corrected' into long format with column for cell line, perturbation, and proliferation, i.e.:
#        Example:
#
#            cell line,perturbation,proliferation
#            T98G.1,miR-34a mimic,1234
#            ...
#
#     b. Make three plots, one for each cell line, with violin plots for each perturbation
#     c. Be sure to include legend, axes labels, and title
#     d. Save as all plots into one PDF called 'perturbation_by_cell_line.pdf'


# Reset index so perturbation becomes a column
prolif_plot_long = prolif_corrected.reset_index().rename(columns={'index': 'perturbation'})

# remove the no cell controls as we don't want them in our plot
prolif_plot_long = prolif_plot_long[~prolif_plot_long['perturbation'].str.startswith('No_Cell')]

replicates = ['HA.1', 'HA.2', 'HA.3', 'HA.4', 'T98G.1', 'T98G.2', 'T98G.3','T98G.4', 'U251.1', 'U251.2', 'U251.3', 'U251.4']

# Melt the replicate columns into long format
prolif_plot_long = pd.melt(prolif_plot_long, id_vars='perturbation', value_vars = replicates, var_name = 'variable', value_name='value')

for index in prolif_plot_long.index:
    row = prolif_plot_long.loc[index]
    if not row['perturbation'].startswith("No_Cell"):
        old_name = prolif_plot_long.loc[index]['variable']
        new_name = old_name[0:-2]
        prolif_plot_long.at[index, 'variable'] = new_name



# Multipage PDF

from matplotlib.backends.backend_pdf import PdfPages

with PdfPages('violinPlots_byCellLine_byPerturbation.pdf') as pdf:
    for cellline in celllines:
        # Filter for current cell line
        subset = prolif_plot_long[prolif_plot_long['variable'] == cellline]

        plt.figure(figsize=(12, 6))
        sns.violinplot(x='perturbation', y='value', data=subset, inner='box')
        plt.title(f'Violin Plot for {cellline}')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        pdf.savefig()
        plt.close()


## 5. (20pts) Compare each treatment versus its corresponding control per cell line and per perturbation
#  Deliverables:
#     a. Set up a dictionary to hold the output of the comparisons called 'perturbation_comparisons'
#         - Hint: Need to capture cell line, perturbation, fold-change, T-statistic, and p-value
#        Example:
#                              HA.FC    HA.stat  HA.pvalue   T98G.FC  T98G.stat   T98G.pvalue   U251.FC  U251.stat   U251.pvalue
#         CEBPD_siRNA       0.670732  -2.104883   0.079936  0.836220  -3.646332  1.075329e-02  0.891289  -4.541457  3.927184e-03
#         miR_29a_Mimic     1.433333   2.861776   0.028735  1.208154   2.326915  5.888886e-02  1.385634   4.243691  5.417672e-03
#                                                                   ...
#
#        * Let pandas decide if scientific notation is needed. In this case p-values for T98G and U251 are in scientifc notation
#        - For the column names use [HA, T98G, U251] connected to [FC, stat, pvalue] using a period, e.g. [HA.FC, HA.stat, HA.pvalue]
#
#     b. Extract the proliferation values for the four replicates for the perturbation of interest
#     c. Extract the proliferation values for the four replicates for the corresponding negative control
#     d. Compute the fold-change as the median of the perturbation proliferation values divided by the
#        median of proliferation for the negative control, and store it into the dictionary
#     e. Compute the Student's t-test p-value using the replicate values for the perturbation versus the
#        replicate values for the negative control (don't use the median value, use four replicate
#        values), and store into the dictionary
#     f. Convert the dictionary into a pandas DataFrame called 'perturbation_comparisons_df'
#     g. Write out the DataFrame to a csv file called 'perturbation_comparisons_df.csv'






## 6. (20pts) Biological interpretation that is to be completed in a separate Word document that will be converted into a PDF for submission
#  Deliverables:
#     HINT: These responses are to be completed in a separate Word document that will be converted into a PDF for submission
#     a. Which many treatments increased or decreased proliferation?
#     b. Which differences are statistically significant?
#     c. What are the differences between the human astrocyte (HA) and glioblastoma lines (T98G and U251)?
#     d. Any trends that could reflect tumor biology?
#     e. Which combinations from Table S19 listed above should be prioritized?
#     c. Be sure to include legend, axes labels, and title
#     d. Save as all plots into one PDF called 'perturbation_by_cell_line.pdf'
