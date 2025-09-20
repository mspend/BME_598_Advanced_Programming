##########################################################
## BME 494/598:  In-class activity 3.2 #1               ##
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
import json
import chemparse as cp


## 1. Load up periodic table of elements data
# Deliverables:
#   a. Open the 'periodic_table_lookup.json' file
#   b. Load periodic table of elements data into variable called 'ptoe_all'
#   c. Print the keys for ptoe_all
#   d. Print the key 'order'
#   e. Print the key 'carbon'
#   f. Print the atomic symbol for 'carbon'
#   g. Print the atomic mass for 'carbon'

with open('periodic_table_lookup.json', 'r') as json_file:
    ptoe_all = json.load(json_file)
    
# print(ptoe_all.keys())
# print(ptoe_all['order'])
# print(ptoe_all['carbon'])
# print(ptoe_all['carbon']['symbol'])
# print(ptoe_all['carbon']['atomic_mass'])



## 2. (5pts) Build a dictionary that is keyed off atomic symobl and has atomic mass as a value
#     Molecules will be given as dictionaries with atomic symbol as key and value being the
#     number of atoms. The ptoe_all is keyed off the full atom name, and we want to have the
#     atomic mass keyed off the atomic symbol. {'C': 12.011, ... }
# Deliverables:
#   a. Initialize a new dictionary 'ptoe'
#   b. Iterate through all 'ptoe_all' atoms (easiest to use ptoe_all['order'] as the iterable)
#   c. Add each atom into 'ptoe' with atomic symbol ('symbol') as the key and atomic mass
#      ('atomic_mass') as the value
#   d. Print the atomic mass of 'C' using 'ptoe', and check if equals 12.011
#   e. Print the atomic mass of 'Fe' using 'ptoe', and check if equals 55.8452
#   f. Print the atomic mass of 'U' using 'ptoe', and check if equals 238.028913

ptoe = {}

# test case
# ptoe[ptoe_all['carbon']['symbol']] = ptoe_all['carbon']['atomic_mass']
# print(ptoe)

for element in ptoe_all['order']:
    ptoe[ptoe_all[element]['symbol']] = ptoe_all[element]['atomic_mass']
    
# print(ptoe)

# print(ptoe['C'])
# print(ptoe['Fe'])
# print(ptoe['U'])


# if ptoe['C'] == 12.011:
#     print('Test 1: Pass')
# else:
#     print('Test 1: Fail')
    
    
# if ptoe['Fe'] == 55.8452:
#     print('Test 2: Pass')
# else:
#     print('Test 2: Fail')


# if ptoe['U'] == 238.028913:
#     print('Test 3: Pass')
# else:
#     print('Test 3: Fail')

## 3. (5pts) Define a function 'molecular_weight' to find the molecular weight for molecules
# Deliverables:
#   a. Write a function 'molecular_weight' that:
#      - Takes two arguments 'molecule' and 'ptoe'
#        a. 'molecule' argument is a dictionary of atoms in the molecule with the value being the number of atoms
#        b. 'ptoe' the dictionary you developed above with atomic symbol as keys and atomic mass as values
#      - Use the atomic weights and the molecular composition to compute the molecular weight
#      - Return the molecular weight
#      - Docstrings to document the function
#   b. Print the results of these tests using the following molecules:
#      - Glycine = C2H5NO2
#      - Glucose = C6H12O6
#      - Palmitic acid = C16H32O2
#      - ATP = C10H16N5O13P3
#      - Dichlorodifluoromethane = CCl2F2
#      - Selenocysteine = C3H7NO2Se
#      - Heme B = C34H32O4N4Fe
#   c. Check molecular weights on the web to make sure your code, function, and tests are working correctly


  
 # Create a test dictionary with the number of atoms per element in the molecule   
# num_atoms_glycine = {'carbon': 2, 'hydrogen': 5, 'nitrogen': 1, 'oxygen': 2}
# print(num_atoms_glycine)

# # figure out how the indexing will work
# # mol_weight = num_atoms_glycine['carbon']*ptoe['C'] + num_atoms_glycine['hydrogen']*ptoe['H']


# mol_weight = 0  
# # iterate over the number of unique elements in the molecule
# for element in num_atoms_glycine.keys():
#     element_symbol = ptoe_all[element]['symbol']
#     # multiply the number of atoms of that element by its atomic weight
#     # ptoe is a dicitonary that stores the atomic weight of each element, with the atomic symbol as the key
#     mol_weight += num_atoms_glycine[element]*ptoe[element_symbol]
# print(mol_weight)
            
def molecular_weight(molecule, ptoe):
    '''
    Calculate the molecular weight of any molecule
    
    Parameters:
        molecule: a dictionary with the atomic symbol of each element it the molecule as the keys and the number of atoms as the values
            Such a dictionary can easily be created with the chemparse package 
        ptoe: a dictionary that stores the atomic weight of each element, with the atomic symbol as the key
            
    Multiply the number of atoms of that element by its atomic weight
    
    Return the molecular weight
    
    '''
    mol_weight = 0
   
    for element in molecule.keys():
        atomic_weight = ptoe[element]
        num_elements_in_molecule = molecule[element]
        mol_weight += num_elements_in_molecule * atomic_weight
        # print(mol_weight)
    return mol_weight
    


# Call the function to test it 

glycine = cp.parse_formula('C2H5NO2')
glucose = cp.parse_formula('C6H12O6')
palmitic_acid = cp.parse_formula('C16H32O2')
ATP = cp.parse_formula('C10H16N5O13P3')
dichlorodifluoromethane = cp.parse_formula('CCl2F2')
selenocysteine = cp.parse_formula('C3H7NO2Se')
heme_b = cp.parse_formula('C34H32O4N4Fe')

# print(molecular_weight(glycine, ptoe)) 
# print(molecular_weight(glucose, ptoe))   
# print(molecular_weight(palmitic_acid, ptoe))  
# print(molecular_weight(ATP, ptoe))  
# print(molecular_weight(dichlorodifluoromethane, ptoe))  
# print(molecular_weight(selenocysteine, ptoe))  
# print(molecular_weight(heme_b, ptoe))  

