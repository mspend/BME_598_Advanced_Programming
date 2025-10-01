import scipy.stats as stats
import pandas as pd
import matplotlib.pyplot as plt

# Load data
file_path = '/Users/maurispendlove/ASU Dropbox/Mauri Spendlove/Classes/BME_598_Programming/Assignments/BME_598_Advanced_Programming/module_6/fev_dat.txt'

# pandas will add an index column automatically
fev_data = pd.read_table(file_path, header=0)

# How many null entries do we have per column?
print(fev_data.info())
# Zero null entries

# How many smokers are there?
# yes = 1, no = 0
print(fev_data['smoke'].value_counts())


print(fev_data['age'].max())
print(fev_data['age'].min())


# T-test comparing sex and FEV
# 0 = female, 1 = male
stats.ttest_ind(fev_data.loc[fev_data.loc[:,'sex']==0, 'FEV'], fev_data.loc[fev_data.loc[:,'sex']==1,'FEV'])

# TtestResult(statistic=-5.441214544130701, pvalue=7.495753953787993e-08, df=652.0)


# T-test comparing age and FEV
stats.ttest_ind(fev_data.loc[fev_data.loc[:,'age']<10, 'FEV'], fev_data.loc[fev_data.loc[:,'age']>=10,'FEV'])

# TtestResult(statistic=-22.41527671223248, pvalue=6.102639320681126e-83, df=652.0)


# T-test comparing smoking status and FEV
# yes = 1, no = 0 
stats.ttest_ind(fev_data.loc[fev_data.loc[:,'smoke']==0, 'FEV'], fev_data.loc[fev_data.loc[:,'smoke']==1,'FEV'])

plt.scatter(fev_data['ht'], fev_data['FEV'])
plt.xlabel("Height")
plt.ylabel("FEV")
plt.title("Height vs FEV")
plt.show()