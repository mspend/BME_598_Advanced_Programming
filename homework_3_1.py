##########################################################
## BME 494/598:  Homework 3.1 Python code               ##
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

## Your name
# Mauri Spendlove

## Import required packages
import chemparse as cp
import pandas as pd


## 1. Install chemparse
#  Deliverables:
#     a. In console window run:
#        pip install chemparse
#     b. Restart console:  Consoles tab -> Restart Kernel
#     c. Import 'chemparse' as 'cp'
#     d. Test chemparse using 'C6H12O6'

glucose_dict = cp.parse_formula('C6H12O6')
print(glucose_dict)



## 2. (10pts) Load up 'inclassactivity_3_2_1.py' data and function
#  Deliverables:
#     a. Import your 'inclassactivity_3_2_1' as 'ica321' so you can use the 'molecular_weight' function
#        - Will have to be the full name of the file without the '.py':
#              import inclassactivity_3_2_1 as ica321
#     b. Test whether the imported 'molecular_weight' function works by printing out the molecular weight of glucose 'C6H12O6'

import inclassactivity_3_2_1 as ica321
periodic_table_dict = ica321.ptoe
print(ica321.molecular_weight(glucose_dict, periodic_table_dict))


## 3. (10pts) Read in molecules from 'molecules.txt', a file that has the columns molecule name and
#     molecular formula. Store into a dictionary with the molecule name as keys and molecular
#     formula as the values.
# Deliverables:
#   a. Open the 'molecules.txt' file
#      - Hint:  might be helpful to look at the file in a text editor first! Delimiter?
#   b. Get rid of header
#   c. Iteratively split apart the contents of each line
#   d. Store into 'molecules' dictionary
#      - keys will be molecule name
#      - values will be molecular formula
#   e. Use a print statement to check everything worked as expected

molecules = {}

with open('molecules.txt', 'r') as molecular_formulas:
    # skip header
    next(molecular_formulas)
    # make each line a key:value pair in a dictionary
    for line in molecular_formulas:
        line = line.strip()
        line = line.split('\t')
        if len(line) != 2:
            print("String splitting didn't work; Line is not the correct size")
        line[0] = (line[0]).lower()
        molecules[line[0]] = line[1]

print(molecules)


## 4. (20pts) Compute molecular weights
# Deliverables:
#   a. Iterate over the molecules
#   b. Compute the molecular_weight of each molecule
#   c. Store into 'molecular_weights' dictionary
#      - keys will be molecule name
#      - values will be molecular weight
#      - You must use Worksheet #9 molecular_weight and chemparse parse_formula functions
#   e. Use a print statement to check everything worked as expected

molecular_weights = {}

for molecule in molecules.keys():
    # print(molecule)
       
    # convert the molecular formula in the molecules dictionary into a dictionary following chemparse format, i.e. {'C': 27.0, 'H': 44.0, 'O': 1.0}
    molecular_formula = molecules[molecule]
    
    num_atoms_per_element = cp.parse_formula(molecular_formula)
    # print(num_atoms_per_element)
    
    # calculate molecular_weight
    molecular_weight = ica321.molecular_weight(num_atoms_per_element, periodic_table_dict)
    
    # populate dictionary with molecular weights
    molecular_weights[molecule] = molecular_weight
    
print(molecular_weights)


## 5. (10pts) Read in DrugBank drug names and IDs from 'DrugBank.csv', a file that has DrugBank drug
#     names and unique computer readable DrugBank IDs. Store into a dictionary 'drugs' with the
#     drug names as keys and drug ID as values.
# Deliverables:
#   a. Open the 'DrugBank.csv' file
#      - Hint:  might be helpful to look at the file in a text editor first! Delimiter?
#   b. Get rid of header
#   c. Iteratively split apart the contents of each line
#   d. Store into 'drugs' dictionary
#      - keys will be drug name
#      - values will be drug ID
#   e. Use a print statement to check everything worked as expected

drugs = {}

with open('DrugBank.csv', 'r') as file:
    # skip header
    next(file)
    
    for line in file:
        line = line.strip()
        line = line.split(',',1)
        
        if len(line) != 2:
            print("String splitting didn't work; Line is not the correct size")
        
        line[1] = (line[1]).lower()
        drugs[line[1]] = line[0]

print(drugs)

## 6. (20pts) Determine which drugs have a molecular weight in ‘molecular_weights’ and store into the
#     dictionary ‘drug_weights’.
# Deliverables:
#   a. Iterate over drug names
#   b. Check if the drug name exists in 'molecules'
#   c. If drug name in 'molecules' then add to 'drug_weights' dictionary
#      - keys will be drug name
#      - values will be the drugs molecular weight
#   d. Use a print statement to check everything worked as expected

drug_weights= {}

# if drug is in the dictionary molecules, which we have molecular weights for, add it to a drug weights dictionary

for drug in drugs.keys():
    if drug in molecules:
        drug_weights[drug] = molecular_weights[drug]

print(drug_weights)
        

## 7. (20pts) Identify drugs with molecular weight within 400 +/- 5 g/mol.
# Deliverables:
#   a. Iterate through drugs with molecular weights ('drug_weights')
#   b. Check if drug molecular weight is within range (400 +/- 5 g/mol)
#   c. If passes these filters then add drug to possible_drugs
#   d. Use a print statement to check everything worked as expected

lower_range = 400 - 5
upper_range = 400 + 5

possible_drugs = []

for drug in drug_weights:
    if drug_weights[drug] > lower_range and drug_weights[drug] < upper_range:
        possible_drugs.append(drug)

print(possible_drugs)



## 8. (10pts) Write out CSV file with all possible drug and target combiniations.
#   a. The file that should have the following header:
#        Drug,DrugBank_ID,MW
#    Description of columns:
#      Drug = Drug name
#      DrugBank_ID = ID from DrugBank in 'drugs'
#      MW = molecular weight from 'drug_weights'
#
# Deliverables:
#   a. Iterate through 'possible_drugs'
#   c. Add information to write out into a variable temporary variable 'tmp1'
#   d. Join information by commas and append to variable 'write_me'
#   e. After iterations, open file 'possible_drugs.csv' to write out
#   f. Write out 'write_me' joined by newlines
#   g. Open in a spreadsheet program to ensure all data is included that was requested

possible_drugs_df = pd.DataFrame()

possible_drugs_df['Drug'] = possible_drugs

drug_ids = []

molecular_weights_filtered = []

for drug in possible_drugs:
    drug_id = drugs[drug]
    drug_ids.append(drug_id)
    
    molecular_weight = molecular_weights[drug]
    molecular_weights_filtered.append(molecular_weight)


possible_drugs_df['DrugBank_ID'] = drug_ids
possible_drugs_df['MW'] = molecular_weights_filtered

print(possible_drugs_df)

possible_drugs_df.to_csv('possible_drugs.csv', index=False)


# Other way, not using pandas:
    
    
# write_me = []

# for drug in possible_drugs:
#     drug_id = drugs[drug]
#     molecular_weight = molecular_weights[drug]
#     tmp1 = [drug, drug_id, molecular_weight]
#     write_me.append(tmp1)

# with open('possible_drugs.csv', 'w') as file:
#     file.write('Drug,DrugBank_ID,MW\n')
#     for drug in write_me:
#         drug[2] = str(drug[2])
#         line = ' '.join(drug)
#         line = line + '\n'
#         file.write(line)
        

