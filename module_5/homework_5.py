import pandas as pd
import matplotlib.pyplot as plt

# What is the  impact of age on a set of genes in 
# bladder urothelial carcinoma, lung adenocarcinoma, and skin cutaneous melanoma?

# Are there sex-specific effects (gender) on overall survival (OS.time)

panCancer_phenos = pd.read_csv('module_5/phenotypes_panCancer.csv')














penguins = sns.load_dataset("penguins")

sns.scatterplot(x="flipper_length_mm", 
                y="body_mass_g", 
                hue="species", 
                data=penguins)
plt.show()