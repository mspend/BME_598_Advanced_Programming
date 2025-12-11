##########################################################
## BME 598:  final_coding_exercise.py                   ##
##  ______     ______     __  __                        ##
## /\  __ \   /\  ___\   /\ \/\ \                       ##
## \ \  __ \  \ \___  \  \ \ \_\ \                      ##
##  \ \_\ \_\  \/\_____\  \ \_____\                     ##
##   \/_/\/_/   \/_____/   \/_____/                     ##
## @Developed by: Plaisier Lab                          ##
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

### Free functions that will be helpful!!!

## Convert  to DNA nucleotides
def convert_DNA_to_emitted_indices(rev_comp_seed_seq):
    emitted_indices = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    return np.array([emitted_indices[i] for i in rev_comp_seed_seq])


## Convert  to DNA nucleotides
def convert_to_states(hidden):
    states = ['NM1', 'PSSM0', 'PSSM1', 'PSSM2', 'PSSM3', 'PSSM4', 'PSSM5', 'PSSM6', 'PSSM7', 'NM2']
    return ([states[i] for i in hidden], ''.join(['|' if not states[i] in ['NM1','NM2'] else '_' for i in hidden]))


###############################
### Develop your code below ###
###############################


## 1. (10pts) Prepare an executive summary that clearly states the goals of the proposed analysis
#             and provides a concise overview of the analytical steps to be undertaken 
#     - The executive summary needs to be written in plain text, not as bullet points.
#     - And it should be no more than one paragraph.
#     - The overview of analytical steps should be written as an outline with bullet points.
#     - Provide more detail than the headings of the Deliverables section.



## 2. (10pts) Load up the packages you will need to run
#     A. Install any packages you haven't already
#     B.Restart console: Consoles tab -> Restart Kernel
#     C. Import all the packages
import pandas as pd
import numpy as np
import hmmlearn as hmm
import random
import json


## 3a. (15pts) Get the PSSM for your analyses
#    A. Set 'sid' variable and get the correct PSSM:
#    - (5pts) Set 'sid' variable correctly to your ASU student ID
#    - (5pts) Select the correct PSSM from the json file

# Change sid to be your Student ID as an integer
sid = <student_id_goes_here_as_int>
random.seed(sid)
pssm_number = random.sample(range(1,2139),1)[0]
with open('pssms.json', 'r') as jsonFile:
    pssms = json.load(jsonFile)
pssm = pssms[str(pssm_number)]

#    B. Create an emission probability matrix for the PSSM as a variable called 'emission_prob'
#      - Rows = [A, C, G, T]
#      - Columns = [NM1, PSSM0, PSSM1, PSSM2, PSSM3, PSSM4, PSSM5, PSSM6, PSSM7, NM2]
#      - NM1: All are equally likely, 0.25
#      - PSSM0 to PSSMn: Based on likelihood of nucleotide in PSSM
#      - NM2: All are equally likely, 0.25
#    - (2pts) 'emission_prob' shape is correct
#    - (3pts) 'emission_prob' values are correct


# values from pssms[str[969]]
emission_prob = np.array([
    [0.25, 0.25, 0.25, 0.25],
    [0.10667944, 0.66006833, 0.13455484, 0.09869739],
    [0.1608945 , 0.7473111 , 0.00979773, 0.08199666],
    [0.01350603, 0.09173628, 0.8475282 , 0.04722949],
    [0.09914967, 0.66028365, 0.13047052, 0.11009616],
    [0.05459905, 0.8040019 , 0.10038415, 0.04101489],
    [0.01034981, 0.77183734, 0.14130188, 0.07651097],
    [0.0927581 , 0.76955479, 0.01959754, 0.11808957],
    [0.15271665, 0.00892138, 0.75811998, 0.08024198],
    [0.25, 0.25, 0.25, 0.25]]
    )


## Values I wrote down from class
# emission_prob = np.array([
#     [0.25, 0.25, 0.25, 0.25],
#     [0.8 , 0.04, 0.05, 0.11],
#     [0.11, 0.85, 0.02, 0.02],
#     [0.75, 0.1 , 0.1 , 0.05],
#     [0.01, 0.01, 0.1 , 0.88],
#     [0.02, 0.1 , 0.02, 0.86],
#     [0.1 , 0.6 , 0.2 , 0.1 ],
#     [0.05, 0.8 , 0.05, 0.1 ],
#     [0.99, 0.  , 0.  , 0.01],
#     [0.25, 0.25, 0.25, 0.25]]
#     )


# pre-transposed:
# emission_prob = np.array([ 
#     [0.25, 0.80, 0.11, 0.75, 0.01, 0.02, 0.10, 0.05, 0.99, 0.25],
#     [0.25, 0.04, 0.85, 0.10, 0.01, 0.10, 0.60, 0.80, 0.00, 0.25],
#     [0.25, 0.05, 0.02, 0.10, 0.10, 0.02, 0.20, 0.05, 0.00, 0.25],
#     [0.25, 0.11, 0.02, 0.05, 0.88, 0.86, 0.10, 0.10, 0.01, 0.25]
# ])




## 4. (10pts) Prepare parameters for a miRvestigator HMM
#    A. Define starting probabilities as the variable 'start_prob'
#    - Order is:   [NM1, PSSM0, PSSM1, PSSM2, PSSM3, PSSM4, PSSM5, PSSM6, PSSM7, NM2]
#      - Equal likelihood to drop in at NM1 or any of the PSSMn states, i.e., 1 out of 9 chances
#      - Zero likelihood to drop into NM2 state
#    - (2pts) start_prob shape
#    - (3pts) start_prob values are correct
#    B. Define transition probabilities as the variable 'trans_prob'
#      - NM1: Staying in NM1 is very unlikely, and transitioning out to the PSSMs is equally likely, but not to NM2
#      - PSSM0 to PSSMn-1: Transition out of PSSMt to PSSMt+1 is highly likely, and to NM2 is unlikely 0.01
#      - PSSMn: Transition to NM2 is 1
#      - NM2: Transition to NM2 is 1
#    - (2pts) trans_prob shape
#    - (3pts) trans_prob values are correct

start_prob = np.array([1/9, 1/9, 1/9, 1/9, 1/9, 1/9, 1/9, 1/9, 1/9, 0])

trans_prob = np.array([
    [0.01, 0.99/8, 0.99/8, 0.99/8, 0.99/8, 0.99/8, 0.99/8, 0.99/8, 0.99/8, 0], #NM1
    [0, 0.99, 0, 0, 0, 0, 0, 0, 0, 0.01], #PSSM0
    [0, 0, 0.99, 0, 0, 0, 0, 0, 0, 0.01], #PSSM1
    [0, 0, 0, 0.99, 0, 0, 0, 0, 0, 0.01], #PSSM2
    [0, 0, 0, 0, 0.99, 0, 0, 0, 0, 0.01], #PSSM3
    [0, 0, 0, 0, 0, 0.99, 0, 0, 0, 0.01], #PSSM4
    [0, 0, 0, 0, 0, 0, 0.99, 0, 0, 0.01], #PSSM5
    [0, 0, 0, 0, 0, 0, 0, 0.99, 0, 0.01], #PSSM6
    [0, 0, 0, 0, 0, 0, 0, 0, 0.99, 0.01], #PSSM7
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1], #NM2
])


## 5. (15pts) Build an miRvestigator HMM instance and set parameters
#    A. Determine the HMM type and instantiate it as 'model1' with the correct number of components and features
#      - Also, set the 'random_state' parameter to 42 
#    - (5pts) Number of n_components
#    - (5pts) Number of n_features
#    B. Set the starting probabilities
#    - (1pts) Check if startprob_ attribute set
#    C. Set the transition probabilities
#    - (2pts) Check if transmat_ attribute set
#    D. Set the emission probabilites
#    - (2pts) Check if emissionprob_ attribute set



## 6. (10pts) Load all the miRNA seed sequences and compare against PSSM miRvestigator model to find best fit
#   A. Load all miRNA seed sequences from the 'miRBase_miRNAs.csv' file into a Pandas DataFrame called 'miRNAs'
#     - (5pts) Check if shape of miRNAs Pandas DataFrame is correct
#   B. Iteratively compare each miRNA seed sequence (or reverse complement, you figure it out) to the miRvestigator HMM
#     - Note:  You may need to reverse complement the seed sequence, also you will need to convert it to state indices
#       - A function is provided to help with this:  'convert_DNA_to_emitted_indices'
#     - The decode function computes the posterior probability and hidden states from running the Viterbi algorithm
#       - Strongly suggest using model1.decode(X), where X is the output of 'convert_DNA_to_emitted_indices'
#       - The decode function has two output variables, store them as 'posterior_prob' and 'hidden':
#           - i.e.,   posterior_prob, hidden = model1.decode(X)
#       - Next, it is strongly suggested that the function 'convert_to_states' is applied to the hidden states, to translate them into something more useful:
#           - i.e.,  states, pattern = convert_to_states(hidden)
#     - Store results from each miRNA comparison into a dictionary and then append it to a list called 'results'
#       - For each comparison, please capture the following as key-value pairs for each dictionary:
#           - 'miRNA' - miRNA name, e.g., hsa-miR-33a-3p
#           - 'seed seq' - the first eight base pairs (bps) of the miRNA that are most involved in determining the RISC complex binding site
#           - 'seed rev comp seq' - the reverse complement of the first eight base pairs (bps) of the miRNA that are most involved in determining the RISC complex binding site
#           - 'states' - the 'states' output of convert_to_states, where the labels of the hidden states for the best fit are described
#           - 'pattern' - the 'pattern' output of convert_to_states, where | = matches, and _ = mismatches
#           - 'matching' - count up the number of matches (|) in pattern
#           - 'posterior_prob' - the log likelihood of the fit for a given miRNA seed sequence versus the miRvestigator model for a PSSM
#       - (5pts) Check length of results object



## 7. (20pts) Write out csv file sorted by descending posterior_prob
#   A. Turn 'results' into a Pandas DataFrame called 'results_df'
#   B. Set 'miRNA' as the index for the DataFrame
#     - (5pts) Check shape of results object
#   C. Sort the DataFrame by 'posterior_prob' so the values are descending from largest to smallest, values are log likelihoods
#     - (10pts) Is top miRNA in 'results_df' the correct one
#   D. Write out the sorted DataFrame as 'results.csv'
#     - (5pts) Check for 'results.csv file'



## 8. (10pts) Biological interpretation that is to be completed in a separate Word document that will be converted into a PDF for submission
#   A. Using a minimum of 3 and a maximum of 5 paragraphs of plain text (not bullet points), describe how the data in the figures and output from statistical tests in the results tables answer the following questions:
#     - Using the results from the comparisons to all miRNAs, which miRNA is the best fit and why?
#     - If the PSSM was identified from the 3' UTRs of a set of co-expressed genes, what does the matching of that miRNA mean for those genes and for their expression?
#   B. You must embed excerpts of the table to reinforce your interpretation.
#     - Tables need to have headers and column headers that are intuitive and descriptive.
#     - Every table must have a callout in the text, e.g., (Table 1)


