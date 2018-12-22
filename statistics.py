from scipy import stats
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests

data = pd.read_csv('cg_target_seq_counts.csv')


def pairwise_t_test(data):

    """ checks anova significance, runs pairwise t test and corrects with Bonferroni's """

    anova = stats.f_oneway(data['target_per_bp'], data['up_target_per_bp'], data['down_target_per_bp'])
    print(anova)

    target_up_t_test = stats.ttest_ind(data['target_per_bp'], data['up_target_per_bp'])
    target_down_t_test = stats.ttest_ind(data['target_per_bp'], data['down_target_per_bp'])
    up_down_t_test = stats.ttest_ind(data['up_target_per_bp'], data['down_target_per_bp'])

    p_values = [target_up_t_test[1], target_down_t_test[1],  up_down_t_test[1]]

    p_adjusted = multipletests(p_values, method='bonferroni')

    print(p_adjusted)

    return p_adjusted

pairwise_t_test(data)