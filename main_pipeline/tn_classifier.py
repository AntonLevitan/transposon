import matplotlib
import matplotlib.pyplot
import numpy
import scipy
import sklearn.ensemble
import sklearn.metrics
import sklearn.model_selection
import pandas as pd
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

# original_training = pd.read_csv('glabrata_curated_labels.csv')
# original_training = pd.read_csv('albicans_training_set.csv')
original_training = pd.read_csv('sp_training.csv')
original_training['label'] = numpy.where(original_training['label'] == 'Essential', 1, 0)

data = pd.read_csv('SP_PB_SRR089408_all_lengths.csv')

data_with_labels = original_training.merge(data, on='id', how='left')
data_with_labels = data_with_labels[data_with_labels['target_seq_counts'] > 10].reset_index()

selected_features = ['Hits', 'Reads', 'Length', 'Neighbourhood Index',# 'Neighbourhood Index 2',
                     'Hits_600_bp_upstream', #'Hits_600_bp_downstream',
                     #'600_bp_up_target_seq_counts', '600_bp_down_target_seq_counts',
                     'Freedom Index', 'target_per_bp']#,	'600_up_target_per_bp',	'600_down_target_per_bp']

features = data_with_labels[selected_features]
labels = data_with_labels['label']

# features = original_training[selected_features]
# labels = original_training['label']
title = 'Sp over 10 target sequence'


def cross_validate(X, y, title):

    # this is where the model is defined and trained, random_state is specified to ensure reproducibility
    cv = sklearn.model_selection.StratifiedKFold(n_splits=5, random_state=0)
    classifier = sklearn.ensemble.RandomForestClassifier(n_estimators=200, n_jobs=-1, random_state=0)
    probas_ = sklearn.model_selection.cross_val_predict(classifier, X, y, method="predict_proba",
                                                        cv=cv, n_jobs=-1)[:, 1]

    # plotting of the ROC taken from the scikit-learn documentation
    tprs = []
    aucs = []
    mean_fpr = numpy.linspace(0, 1, 100)

    i = 0
    for train, test in cv.split(X.values, y.values):
        fpr, tpr, thresholds = sklearn.metrics.roc_curve(y[test], probas_[test])
        tprs.append(scipy.interp(mean_fpr, fpr, tpr))
        tprs[-1][0] = 0.0
        roc_auc = sklearn.metrics.auc(fpr, tpr)
        aucs.append(roc_auc)
        matplotlib.pyplot.plot(fpr, tpr, lw=1, alpha=0.3, label="ROC fold %d (AUC = %0.2f)" % (i, roc_auc))

        i += 1
    pd.DataFrame([fpr, tpr, thresholds]).T.to_csv(title + 'tprs.csv')
    matplotlib.pyplot.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r', label="Luck", alpha=.8)

    mean_tpr = numpy.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = sklearn.metrics.auc(mean_fpr, mean_tpr)
    std_auc = numpy.std(aucs)
    matplotlib.pyplot.plot(mean_fpr, mean_tpr, color='b', label=r"Mean ROC (AUC = %0.3f $\pm$ %0.3f)" %
                                                                (mean_auc, std_auc), lw=2, alpha=.8)

    std_tpr = numpy.std(tprs, axis=0)
    tprs_upper = numpy.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = numpy.maximum(mean_tpr - std_tpr, 0)
    matplotlib.pyplot.fill_between(mean_fpr, tprs_lower, tprs_upper, color="grey", alpha=.2, label=r"$\pm$ 1 std. dev.")

    matplotlib.pyplot.xlim([-0.05, 1.05])
    matplotlib.pyplot.ylim([-0.05, 1.05])
    matplotlib.pyplot.xlabel("False Positive Rate")
    matplotlib.pyplot.ylabel("True Positive Rate")
    matplotlib.pyplot.title("Receiver operating characteristic")
    matplotlib.pyplot.legend(loc="lower right")
    matplotlib.pyplot.savefig(title + '.png')
    matplotlib.pyplot.show()

    data_with_labels['probability'] = probas_
    data_with_labels.to_csv(title + '.csv')
    # original_training['probability'] = probas_
    # original_training.to_csv(title + '.csv')

    classifier.fit(X, y)

    feature_list = selected_features
    # feature_list = list(features.columns.values)
    # feature_list = ['Hits', 'Reads', 'Length', 'Neighbourhood Index',# 'Neighbourhood Index 2',
    #                  'Hits_600_bp_upstream',# 'Hits_600_bp_downstream',
    #                  '600_bp_up_target_seq_counts', '600_bp_down_target_seq_counts',
                     # 'Freedom Index', 'target_per_bp']#,	'600_up_target_per_bp',	'600_down_target_per_bp']

    importances = list(classifier.feature_importances_)
    feature_importances = [(feature, round(importance, 2)) for feature,
                                                               importance in zip(feature_list, importances)]
    feature_importances = sorted(feature_importances, key=lambda x: x[1], reverse=True)
    pd.DataFrame(feature_importances).to_csv(title + ' feature importances.csv')

    matplotlib.pyplot.style.use('fivethirtyeight')
    x_values = list(range(len(importances)))
    matplotlib.pyplot.bar(x_values, importances, orientation='vertical')
    matplotlib.pyplot.xticks(x_values, feature_list, rotation=90)
    matplotlib.pyplot.ylabel('Importance')
    matplotlib.pyplot.title('Feature Importance')
    matplotlib.pyplot.savefig(title + '_feature_importance' + '.png')
    matplotlib.pyplot.show()


cross_validate(features, labels, title)
