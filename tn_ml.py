import matplotlib
import matplotlib.pyplot
import numpy
import scipy
import sklearn.ensemble
import sklearn.metrics
import sklearn.model_selection
import pandas as pd

CV_FILE_SUFFIX = ".png"

data = pd.read_csv('/home/user/Desktop/transposon/FreadI1_trimmed_sorted_analysis.csv', header=0)

selected_features = [1, 3, 4, 5, 6, 7, 9]
features = data.iloc[:, selected_features]

ground_truth = pd.read_csv('training_set_Sc.csv', header=0)
ground_truth = ground_truth.iloc[:, [0, 3]]
ground_truth = ground_truth.replace('Essential', 1)
ground_truth = ground_truth.replace('Non essential', 0)

final_data = ground_truth.merge(features)

label = final_data.iloc[:, 1]
final_features = final_data.iloc[:, [2, 3, 4, 5, 6, 7]]

def cross_validate(X, y):
    
    # this is where the model is defined and trained, random_state is specified to ensure reproducibility
    cv = sklearn.model_selection.StratifiedKFold(n_splits=5, random_state=15)
    classifier = sklearn.ensemble.RandomForestClassifier(n_jobs=-1, random_state=15)
    probas_ = sklearn.model_selection.cross_val_predict(classifier, X, y, method="predict_proba", cv=cv, n_jobs=6)[:, 1]

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
    matplotlib.pyplot.plot(mean_fpr, mean_tpr, color='b', label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc), lw=2, alpha=.8)

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


cross_validate(final_features, label)
