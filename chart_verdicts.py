import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr

data = pd.read_excel('sanity_check_2.xlsx')

acds = data.iloc[:, [0, 2]]
hermes = data.iloc[:, [4, 5]]
hermes.columns = ['Standard name', 'Hermes Verdict']
merged = acds.merge(hermes)
# print(merged)
merged.to_csv('merged.csv')

print(pearsonr(merged['AcDs Verdict'], merged['Hermes Verdict']))

sanity_check = pd.read_excel('sanity_check.xlsx').dropna()
print(pearsonr(sanity_check['Sc RF verdict ScTn'], sanity_check['Sp RF verdict SpTn']))

# plt.plot(sanity_check['Sc RF verdict ScTn'], sanity_check['Sp RF verdict SpTn'], 'o', color='green')
# plt.savefig('test14.png')

plt.plot(list(merged['AcDs Verdict']), list(merged['Hermes Verdict']), 'o', color='black', ms=2)
plt.xlabel('AcDs')
plt.ylabel('Hermes')
plt.savefig('test19.png')
plt.show()