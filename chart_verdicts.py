import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr
import numpy as np
from matplotlib.ticker import NullFormatter

data = pd.read_excel('sanity_check_6.xlsx')

acds = data.iloc[:, [0, 2]]
hermes = data.iloc[:, [4, 5]]
hermes.columns = ['Standard name', 'Hermes Verdict']
merged = acds.merge(hermes)
data = merged

# data = pd.read_excel('/home/user/Desktop/transposon/verdict_difference.xlsx', header=0)
# data = data.sort_values(by=['verdict difference'])

# data = pd.read_excel('swathi data.xlsx', header=0)
print(data.head())

tested_feature_1 = 'AcDs Verdict'
tested_feature_2 = 'Hermes Verdict'

print(pearsonr(data[tested_feature_1], data[tested_feature_2]))

x = data[tested_feature_1]
y = data[tested_feature_2]

# definitions for the axes
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
bottom_h = left_h = left + width + 0.02

rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.2]
rect_histy = [left_h, bottom, 0.2, height]

# start with a rectangular Figure
plt.figure(1, figsize=(8, 8))
axScatter = plt.axes(rect_scatter)
plt.xlabel(tested_feature_1, fontsize=18)
plt.ylabel(tested_feature_2, fontsize=18)

axHistx = plt.axes(rect_histx)
axHisty = plt.axes(rect_histy)

# no labels
nullfmt = NullFormatter()
axHistx.xaxis.set_major_formatter(nullfmt)
axHisty.yaxis.set_major_formatter(nullfmt)

# the scatter plot:
axScatter.scatter(x, y, s=6)

xbins = np.linspace(min(x), max(x), 30)
ybins = np.linspace(min(y), max(y), 30)

axHistx.hist(x, bins=xbins)
axHisty.hist(y, bins=ybins, orientation='horizontal')

plt.savefig('both_' + tested_feature_1 + ' vs. ' + tested_feature_2 + '.png')
plt.show()