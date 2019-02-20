import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('filtering_target_seqs_PB_Ca.csv')
x = data['target sequences cutoff']
genes = data['genes']
auc = data['auc']
target = data['Target / bp FI']

fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('Number of target sequences in ORF used as cutoff')
ax1.set_ylabel('AUC', color=color)
ax1.scatter(x, auc, color=color)
for i, txt in enumerate(genes):
    ax1.annotate(txt, (x[i], auc[i]))
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Targets / bp Feature Importance', color=color)  # we already handled the x-label with ax1
ax2.scatter(x, target, color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()
plt.savefig('Ca PB AUC with different targets cutoffs.png')  # otherwise the right y-label is slightly clipped
plt.show()
