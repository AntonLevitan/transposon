import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr

data = pd.read_excel('sanity_check_2.xlsx')

acds = data.iloc[:, [0, 2]]
hermes = data.iloc[:, [0, 1]]
hermes.columns = ['Standard name', 'Anton Verdict']
merged = acds.merge(hermes)
# print(merged)
merged.to_csv('merged.csv')

print(pearsonr(merged['Vlad Verdict'], merged['Anton Verdict']))

sanity_check = pd.read_excel('sanity_check.xlsx').dropna()
print(pearsonr(sanity_check['Sc RF verdict ScTn'], sanity_check['Sp RF verdict SpTn']))

# plt.plot(sanity_check['Sc RF verdict ScTn'], sanity_check['Sp RF verdict SpTn'], 'o', color='green')
# plt.savefig('test14.png')

plt.plot(list(merged['Vlad Verdict']), list(merged['Anton Verdict']), 'o', color='green')
plt.savefig('test18.png')
plt.show()