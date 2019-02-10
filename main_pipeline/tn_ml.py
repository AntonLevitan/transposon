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

CV_FILE_SUFFIX = ".png"

data2 = pd.read_csv('all_genes_sc_ml.csv')
data = pd.read_csv('Both_ml.csv')
# data = pd.read_excel('/home/user/Desktop/transposon/sanity_check_classifier.xlsx', header=0)

selected_features = ['Standard name', 'id', 'Hits', 'Reads', 'Length', 'Neighbourhood Index',
                     '600_bp_up_target_seq_counts', '600_bp_down_target_seq_counts',
                     'Freedom Index', 'target_per_bp',	'600_up_target_per_bp',	'600_down_target_per_bp']

features = data[selected_features]
all_genes_features = features[selected_features[2:]]
features2 = data2[selected_features]
all_genes_features2 = features2[selected_features[2:]]

ground_truth = pd.read_csv('training_set_Sc.csv', header=0)
ground_truth = ground_truth.iloc[:, [0, 3]]
ground_truth = ground_truth.replace('Essential', 1)
ground_truth = ground_truth.replace('Non essential', 0)


glab_train = pd.read_csv('curation.csv', header=0)


final_data = ground_truth.merge(features)
final_data2 = glab_train.merge(features2)

label = final_data['Final usage in training']
final_features = final_data[selected_features[2:]]

label2 = glab_train['label']
final_features2 = glab_train[selected_features[2:]]


def cross_validate(X, y):
    
    # this is where the model is defined and trained, random_state is specified to ensure reproducibility
    cv = sklearn.model_selection.StratifiedKFold(n_splits=5, random_state=0)
    classifier = sklearn.ensemble.RandomForestClassifier(n_estimators=200, n_jobs=-1, random_state=0)
    probas_ = sklearn.model_selection.cross_val_predict(classifier, X, y, method="predict_proba",
                                                        cv=cv, n_jobs=-1)[:, 1]
    classifier.fit(X, y,)
    x = classifier.predict_proba(all_genes_features2)[:, 1]
    # glab_train['Hermes probability'] = x
    # glab_train.to_csv('glabrata_training_verdict.csv')
    prob = pd.Series(x)
    # verdict = pd.DataFrame()
    # verdict['id'] = data2['id']
    # # verdict['label'] = y
    # verdict['Hermes probability'] = prob
    # verdict.to_csv('glabrata_all_verdict.csv')
    data2['Hermes probability'] = prob
    data2.to_csv('sc_all_data.csv')
    metrics = pd.DataFrame()
    fpr, tpr, thresholds = sklearn.metrics.roc_curve(y, probas_)
    metrics['FPR'] = list(fpr)
    metrics['TPR'] = list(tpr)
    metrics['Thresholds'] = list(thresholds)
    # metrics.to_csv('ml_metrics_glabrata.csv')
    # plotting of the ROC taken from the scikit-learn documentation
    tprs = []
    aucs = []
    mean_fpr = numpy.linspace(0, 1, 100)

    i = 0
    for train, test in cv.split(X, y):
        fpr, tpr, thresholds = sklearn.metrics.roc_curve(y[test], probas_[test])
        tprs.append(scipy.interp(mean_fpr, fpr, tpr))
        tprs[-1][0] = 0.0
        roc_auc = sklearn.metrics.auc(fpr, tpr)
        aucs.append(roc_auc)
        matplotlib.pyplot.plot(fpr, tpr, lw=1, alpha=0.3, label="ROC fold %d (AUC = %0.2f)" % (i, roc_auc))

        i += 1

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
    matplotlib.pyplot.savefig("glabrata" + CV_FILE_SUFFIX)
    matplotlib.pyplot.show()


def plot_feat_imp():

    # feature importance plot
    classifier = sklearn.ensemble.RandomForestClassifier(n_estimators=200, n_jobs=-1, random_state=15)
    classifier.fit(final_features2, label2)
    
    # feature_list = list(final_features.columns.values)
    feature_list = ['Hits', 'Reads', 'Length', 'Neighborhood Index',
                     '600_bp_up_hits', '600_bp_down_hits',
                     'Freedom Index', 'target_per_bp',	'600_up_targets',	'600_down_targets']

    importances = list(classifier.feature_importances_)
    feature_importances = [(feature, round(importance, 2)) for feature, 
                           importance in zip(feature_list, importances)]
    feature_importances = sorted(feature_importances, key=lambda x: x[1], reverse=True)
    
    matplotlib.pyplot.style.use('fivethirtyeight')
    x_values = list(range(len(importances)))
    matplotlib.pyplot.bar(x_values, importances, orientation='vertical')
    matplotlib.pyplot.xticks(x_values, feature_list, rotation=90)
    matplotlib.pyplot.ylabel('Importance')
    matplotlib.pyplot.title('Feature Importance')
    matplotlib.pyplot.savefig('glabrata' + '_feature_importance')
    matplotlib.pyplot.show()


def cross_species_validate(X, y):
    classifier = sklearn.ensemble.RandomForestClassifier(n_estimators=200, n_jobs=-1, random_state=0)
    classifier.fit(X, y,)
    probas_ = classifier.predict_proba(final_features2)[:, 1]

    tprs = []
    aucs = []
    mean_fpr = numpy.linspace(0, 1, 100)

    fpr, tpr, thresholds = sklearn.metrics.roc_curve(label2, probas_)
    tprs.append(scipy.interp(mean_fpr, fpr, tpr))
    tprs[-1][0] = 0.0
    roc_auc = sklearn.metrics.auc(fpr, tpr)
    aucs.append(roc_auc)

    matplotlib.pyplot.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r', label="Luck", alpha=.8)

    mean_tpr = numpy.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = sklearn.metrics.auc(mean_fpr, mean_tpr)
    matplotlib.pyplot.plot(mean_fpr, mean_tpr, color='b', label=r"ROC (AUC = %0.3f)" % mean_auc, lw=2, alpha=.8)
    matplotlib.pyplot.xlim([-0.05, 1.05])
    matplotlib.pyplot.ylim([-0.05, 1.05])
    matplotlib.pyplot.xlabel("False Positive Rate")
    matplotlib.pyplot.ylabel("True Positive Rate")
    matplotlib.pyplot.title("Receiver operating characteristic")
    matplotlib.pyplot.legend(loc="lower right")
    matplotlib.pyplot.savefig("sc on glabrata" + CV_FILE_SUFFIX)
    matplotlib.pyplot.show()


# plot_feat_imp()
cross_validate(final_features, label)
# cross_species_validate(final_features, label)
