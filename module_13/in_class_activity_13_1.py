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
# For PHX, AZ, it is very likely to start out sunny. 
# PHX has ~34 days of rain per year
start = { "sun": 0.95, "rain": 0.05 }


## 3. Define transition probabilities
# These transitions mean that if it's sunny, it's very likely to stay sunny.
# If it's rainy, its slightly more likely to switch back to sunny. 
# Often in Phoenix, we only have 1 day of rain at a time. We don't often have multiple days of rain on end.
transitions = {"sun": { "sun": 0.95, "rain": 0.05},
               "rain": {"sun": 0.55, "rain": 0.45}}


## 4. MarkovChain function to predict states
def MarkovChain(start, transitions, num_states):
    # State chain
    states = []
        
    # Ensure the probabilities sum to 1 (a good practice check)
    if not np.isclose(sum(start.values()), 1.0):
        raise ValueError("Start likelihoods must sum to approximately 1.0")
    
    # First, choose a starting state 
    # use the initial probabilities to randomly choose a starting state
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


# ## 6. Sample 50 states from chain
# print(result.sample(50))

# Result:
# I ran this model half a dozen times and got the following numbers of days of rain:
# 28, 35, 26, 22, 22, 25, 47, 43, 44, 35, 40, 30

num_rainy_days = np.array([28, 35, 26, 22, 22, 25, 47, 43, 44, 35, 40, 30])
mean_value = np.mean(num_rainy_days)
print(f"The average number of days of rain: {mean_value}")