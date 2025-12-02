##########################################################
## OncoMerge:  in_class_activity_13_1.py                ##
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


## 1. Load libraries
import numpy as np


## 2. Define starting probabilities
# In our weather example, the 2 possible states are sunny and rainy.
# We are starting with the assumption that these two states are equally likely.
start = { "sun": 0.5, "rain": 0.5 }


## 3. Define transition probabilities
# These transitions mean that if you're already in one weather state, it's more likely to change that weather 
transitions = {"sun": { "sun": 0.8, "rain": 0.2},
               "rain": {"sun": 0.3, "rain": 0.7}}


## 4. MarkovChain function to predict states
def MarkovChain(start, transitions, num_states):
    # State chain
    states = []
        
    # Ensure the probabilities sum to 1 (a good practice check)
    if not np.isclose(sum(start.values()), 1.0):
        raise ValueError("Start likelihoods must sum to approximately 1.0")
    
    # First, choose a starting state using 
    # uses the initial probabilities to randomly choose a starting state
    states.append(np.random.choice(list(start.keys()), size=1, p=list(start.values()))[0])
    
    # Next add other states
    for i in range(num_states-1):
        states.append(np.random.choice(list(transitions[states[-1]].keys()), size=1, p=list(transitions[states[-1]].values()))[0])
    
    # Return states
    return(states)


## 5. Compute 50 samples
# We predict states for an entire year, 365 days
# # of states needs to be an integer. Can'tdo 365.25 states.
result = MarkovChain(start, transitions, 365)
print(result)
print('sun = '+str(np.sum([i=='sun' for i in result])))
print('rain = '+str(np.sum([i=='rain' for i in result])))


## 6. Sample 50 states from chain
print(model.sample(50))

