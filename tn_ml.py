import matplotlib
import matplotlib.pyplot
import numpy
import scipy
import sklearn.ensemble
import sklearn.metrics
import sklearn.model_selection
import pandas as pd
from sklearn import preprocessing
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

CV_FILE_SUFFIX = ".png"

# data = pd.read_csv('/home/user/Desktop/transposon/FreadI1_trimmed_sorted_analysis.csv', header=0)
data = pd.read_excel('/home/user/Desktop/transposon/sanity_check_classifier.xlsx', header=0)

# selected_features = [1,3,4,5,6,7,8]
selected_features = [0,1,2,3,4,5,6]

features = data.iloc[:, selected_features]
all_genes_features = features.iloc[:, [1,2,3,4,5,6]]
# print(all_genes_features)

ground_truth = pd.read_csv('training_set_Sc.csv', header=0)
ground_truth = ground_truth.iloc[:, [0, 3]]
ground_truth = ground_truth.replace('Essential', 1)
ground_truth = ground_truth.replace('Non essential', 0)

final_data = ground_truth.merge(features)

label = final_data.iloc[:, 1]
final_features = final_data.iloc[:, [2, 3, 4, 5, 6, 7]]
# final_features = preprocessing.normalize(final_features)
# print(final_features)

def cross_validate(X, y, data, all_genes_features, all_genes):
    
    # this is where the model is defined and trained, random_state is specified to ensure reproducibility
    cv = sklearn.model_selection.StratifiedKFold(n_splits=5, random_state=0)
    classifier = sklearn.ensemble.RandomForestClassifier(n_estimators=100, n_jobs=-1, random_state=0)
    probas_ = sklearn.model_selection.cross_val_predict(classifier, X, y, method="predict_proba", cv=cv, n_jobs=-1)[:, 1]
    # print(list(probas_))
    classifier.fit(X, y,)
    x = classifier.predict_proba(all_genes_features)[:, 1]
    # print(all_genes_features)
    # print(len(x))
    all_genes['Hermes Verdict'] = x
    all_genes.to_csv('all_anton_verdict.csv')
    # print(all_genes)
    prob = pd.Series(probas_)
    verdict = pd.DataFrame()
    verdict['name'] = data['Standard name']
    verdict['label'] = y
    verdict['verdict'] = prob
    verdict.to_csv('training_hermes_verdict.csv')
    
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
    matplotlib.pyplot.plot(mean_fpr, mean_tpr, color='b', label=r"Mean ROC (AUC = %0.3f $\pm$ %0.3f)" % (mean_auc, std_auc), lw=2, alpha=.8)

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
    matplotlib.pyplot.savefig("sc_tn" + CV_FILE_SUFFIX)
    matplotlib.pyplot.show()


cross_validate(final_features, label, ground_truth, all_genes_features, data)

def plot_feat_imp():
    #feature importance plot
    classifier = sklearn.ensemble.RandomForestClassifier(n_estimators=200, n_jobs=-1, random_state=15)
    classifier.fit(final_features, label)
    
    feature_list=list(final_features.columns.values)
    importances = list(classifier.feature_importances_)
    feature_importances = [(feature, round(importance, 2)) for feature, 
                           importance in zip(feature_list, importances)]
    feature_importances = sorted(feature_importances, key = lambda x: x[1], reverse = True)
    
    # fig,ax = matplotlib.pyplot.subplots(figsize=(8, 5))

    matplotlib.pyplot.style.use('fivethirtyeight')
    x_values = list(range(len(importances)))
    matplotlib.pyplot.bar(x_values, importances, orientation = 'vertical')
    matplotlib.pyplot.xticks(x_values, feature_list, rotation=45)
    matplotlib.pyplot.ylabel('Importance'); matplotlib.pyplot.title('Feature Importance');
    #matplotlib.pyplot.gcf().sublots_adjust(bottom=0.15)
    matplotlib.pyplot.savefig('sc_tn' + '_feature_importance')
    matplotlib.pyplot.show()

# plot_feat_imp()