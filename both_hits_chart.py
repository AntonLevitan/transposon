import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr

data = pd.read_excel('both_hits.xlsx').dropna()
# print(data)
# print(pearsonr(data['Hermes Hits'], data['AcDs Hits']))

# plt.plot(data['Hermes Hits'], data['AcDs Hits'], 'o', color='green')
# plt.plot(data['Hermes NI'], data['AcDs NI'], 'o', color='blue')
plt.plot(data['Hermes FI'], data['AcDs FI'], 'o', color='blue')
plt.savefig('test13.png')