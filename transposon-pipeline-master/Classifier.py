import Shared
import GenomicFeatures
import SummaryTable
import Organisms

import glob
import os
import csv
import cPickle
from itertools import product, chain, repeat
import re
from collections import OrderedDict
from pprint import pprint
import shutil

from sklearn.cross_validation import StratifiedKFold
from sklearn.metrics import roc_curve, auc
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier

import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3

def test_classifier(classifier, features, annotations):
    """Test a classifier using given features and annotations. Uses a 5-fold
    cross-validation method to score all of the test data points.
    
    Parameters
    ----------
    classifier : sklearn Classifier
        The classifier to test.
    features : list of sequences
        A list of features per gene.
    annotations : list of 1 or 0
        A list of annotations for the genes.
    
    Returns
    -------
    tuple of tuples
        The ROC curve parameters, as returned by `sklearn.metrics.roc_curve`,
        and the scores, as the second item of the tuple.
    """
    
    scores = np.empty(len(annotations))
    skf = StratifiedKFold(y=annotations, n_folds=5, shuffle=True, random_state=0)
    for roc_train, roc_test in skf:
        roc_train_set = [features[i] for i in roc_train]
        roc_test_set = [features[i] for i in roc_test]
        roc_train_essentiality = [annotations[i] for i in roc_train]
        
        classifier.fit(roc_train_set, roc_train_essentiality)
        scores[roc_test] = classifier.predict_proba(roc_test_set)[:,1]
    
    return roc_curve(annotations, scores), scores

def get_training_groups():
    """Get all interesting permutations of features to train on.
    
    Note
    ----
    
    Use with caution, as this can be quite large.
    
    Returns
    -------
    generator of tuples
        Each tuple contains a permutation of features to use for training.
    """
    
    subgroups = (("neighborhood_index", "insertion_index"),
                 ("length", "hits", ("length", "hits",), ""),
                 ("max_free_region", "freedom_index", "adj_max_free_region", "adj_freedom_index", ""),
                 ("upstream_hits_50", "upstream_hits_100", ("upstream_hits_50", "upstream_hits_100"), ""),
                 ("std_dev", "index_of_dispersion", ""))
    
    def _flatten(features):
        # Recursive, generator-based flattening.
        # NB: the above is probably a lie, and it's not recursive, but goes
        # one level deep, which is sufficient here.
        for f in features:
            if isinstance(f, str):
                yield f
            else:
                for sub_f in f:
                    yield sub_f
    
    for group in product(*subgroups):
        yield tuple(f for f in _flatten(group) if f)


def get_youden_statistic(fpr, tpr, thresholds):
    """The Youden index (or J statistic) is the threshold which maximizes 
    the (TPR-FPR) value."""
    
    return max(zip(thresholds, tpr, fpr), key=lambda (_th, t, f): t-f)

def draw_venn(title, output_file_path, data, labels):
    """A convenience function for generating Venn diagrams and saving them to
    a destination on the drive. The `data` and `labels` can be of size 2 or 3.
    
    Parameters
    ----------
    title : str
        The title of the plot.
    output_file_path : str
        The path in which the plot will be saved.
    data : list of set
        The data to plot - should be represented as sets.
    labels : list of str
        The labels of the data to plot.
    """
    
    # TODO: should go into a drawing module.
    plt.figure(figsize=(8,8))
    
    if len(data) == 2:
        venn = venn2
    else:
        venn = venn3
    v = venn(data, labels)
    
    for l in v.set_labels:
        l.set_style('italic')
    
    plt.title(title)
    plt.savefig(output_file_path,
                transparent=True,
                dpi=300)
    
    plt.close()
    
def write_group_legends(feature_groups, output_file):
    with open(output_file, 'w') as out_file:
        writer = csv.writer(out_file, delimiter=',')
        writer.writerow(["Group name", "Features"])
        for group_name, columns_for_classification in feature_groups:
            writer.writerow([group_name, ", ".join(columns_for_classification)])

def only_orfs(records):
    return (r for r in records if "ORF" in r["feature"].type)

def no_dubious_features(records):
    return (r for r in records if r["feature"].feature_qualifier != "Dubious")

def no_depleted_hit_features(records):
    # TODO: the deleted genes are known for the cerevisiae (and pombe?) strains, should be added
    # from the list instead of inferred.
    # If there were no insertions in the feature AND its neighborhood,
    # we can't trust this feature (maybe it was deleted):
    return (r for r in records if r["neighborhood_hits"] > 0 or r["hits"] > 0)

def combine_pre_post(pre_records, post_records):
    """Combine pre- and post- records into a single record list, by appending
    prefixes to the record values."""
    
    pre_records_sorted = sorted(pre_records, key=lambda r: r["feature"].standard_name)
    post_records_sorted = sorted(post_records, key=lambda r: r["feature"].standard_name)
    assert len(pre_records_sorted) == len(post_records_sorted)
    
    result = []
    for pre_record, post_record in zip(pre_records_sorted, post_records_sorted):
        assert pre_record["feature"] == post_record["feature"]
        new_record = {}
        for key, value in pre_record.items():
            new_record[key] = new_record["pre_" + key] = value
        for key, value in post_record.items():
            new_record["post_" + key] = value
        result.append(new_record)
            
    return result

def explode_col_config(col_config):
    """"Explodes" (in PHP parlance) the column configuration to display pre-
    and post- data, by appending prefixes to the column names."""
    
    result = []
    
    for config in col_config:
        if config.get("local", False):
            for prefix in ("pre", "post"):
                new_config = dict(config)
                new_config["field_name"] = "%s_%s" % (prefix, new_config["field_name"])
                new_config["csv_name"] = "%s - %s" % (prefix.title(), new_config["csv_name"])
                result.append(new_config)
        else:
            result.append(config)
    
    return result

def run_pipeline(cls_factory, cls_feature_groups, data, train_col_config, class_col_config, output_folder, target_fpr):
    """Classify the given datasets and store the output into the given folder.
    
    Parameters
    ----------
    cls_factory : dict of str to factory function
        A dict of classifier names to their factory functions. Each dataset
        will be classified separately with each of the classifiers given.
        
        For example:
        {
            "Cls1": lambda: MyClassifier(p1=42, p2=666)
        }
    cls_feature_groups : dict of str to tuple
        The feature groups which to use during classification. Each classifier
        will be run with each of the feature groups separately.
        
        For example:
        {
            "Group 1": ("hits", "length", "reads")
        }
    data : dict of str to dataset descriptor
        The datasets to classify. Supports multiple separate classifications.
        Organized around training sets - each training set can be used to
        classify multiple datasets, and be comapred with multiple benchrmarks
        (if desired).
    
        IMPORTANT: records can be mutated (i.e., keys added) during the classification.
        It is recommended to shallow-copy records passed into this function,
        and keep a handle to them. That shallow copy can be kept around for
        further inspection of classifier scores later.
        
        Detailed format:
        {
            "Classificaiton name 1": (
                list_of_essential_records_for_training,
                list_of_non_essential_records_for_training,
                {
                    "Dataset to classify 1": records_to_classify,
                    ...
                },
                # OPTIONAL: dict of benchmarks for comparison
                {
                    "Benchmark 1": (
                        essential_records,
                        non_essential_records
                    ),
                    ...
                }
            ),
            ...
        }
    train_col_config : dict
        TBD:
    class_col_config : dict
        TBD
    output_folder : str
        TBD
    target_fpr : float
        TBD
    
    Returns
    -------
    tuple of tuples
        The ROC curve parameters, as returned by `sklearn.metrics.roc_curve`,
        and the scores, as the second item of the tuple.
    """
    
    # Insert the classifier data into the column configuration.
    # We assume that the last two columns are type and description, so we want them to be last.
    train_col_config = list(train_col_config)
    class_col_config = list(class_col_config)
    
    train_col_config[-2:-2] = \
        [{"field_name": "ground_truth", "csv_name": "Ess. ground truth"}] + \
        sum([[{"field_name": "train-%s-%s" % key, "csv_name": "%s - %s" % key, "format": "%.3f", "local": True},
              {"field_name": "train-%s-%s-verdict" % key, "csv_name": "%s - %s - ess. for FPR %.3f" % (key + (target_fpr,)), "local": True}]
             for key in product(cls_factory.keys(), cls_feature_groups.keys())],
            [])
    
    # If the training set is taken from within the same organism as the unknown
    # genes, we should clearly mark those genes that were used in training. 
    class_col_config[-2:-2] = \
        [{"field_name": "ground_truth", "csv_name": "Ess. ground truth"}] + \
        sum([[{"field_name": "%s-%s" % key, "csv_name": "%s - %s" % key, "format": "%.3f", "local": True},
              {"field_name": "%s-%s-verdict" % key, "csv_name": "%s - %s - ess. for FPR %.3f" % (key + (target_fpr,)), "local": True}]
             for key in product(cls_factory.keys(), cls_feature_groups.keys())],
            [])
    threshold_ixs = {}
    for data_key in data.keys():
        data_folder = os.path.join(output_folder, "classification - %s" % data_key)
        Shared.make_dir(data_folder)
        
        train_ess_records, train_non_ess_records, class_datasets, benchmarks = \
                (data[data_key] + ({},))[:4]
        train_all_records = train_ess_records + train_non_ess_records
        train_annotations = [1] * len(train_ess_records) + [0] * len(train_non_ess_records)
        
        for grp_name, features in cls_feature_groups.items():
            train_all_features = [[r[f] for f in features] for r in train_all_records]
            
            for cls_name, cls_creator in cls_factory.items():
                record_label = "%s-%s" % (cls_name, grp_name)
                train_record_label = "train-" + record_label
                verdict_record_label = record_label + "-verdict"
                train_verdict_record_label = "train-" + record_label + "-verdict"
                
                # Test classifier on training data
                ((fprs, tprs, thresholds), scores) = test_classifier(
                    cls_creator(),
                    train_all_features,
                    train_annotations
                )
                
                # Write out AUC curve
                write_auc_curve(
                    (fprs, tprs, thresholds),
                    os.path.join(data_folder, "train_AUC.%s.png" % record_label)
                )
                
                training_threshold_ix = np.absolute(fprs - target_fpr).argmin()
                threshold_ixs[(cls_name, grp_name)] = training_threshold_ix
                fpr_threshold = thresholds[training_threshold_ix]
                
                for r, s in zip(train_all_records, scores):
                    r[train_record_label] = s
                    r[train_verdict_record_label] = "Yes" if s >= fpr_threshold else "No"

                # Train classifier
                classifier = cls_creator()
                classifier.fit(train_all_features, train_annotations)
                
                # Feature importances?
                if hasattr(classifier, "feature_importances_"):
                    print data_key, cls_name, grp_name
                    pprint(dict(zip(features, classifier.feature_importances_)))
                    print "\n"
                elif hasattr(classifier, "coef_"):
                    print data_key, cls_name, grp_name
                    pprint(dict(zip(features, classifier.coef_[0])))
                    print "\n"
                
                # Classify records for classification, store in records
                for classification_records in chain(*class_datasets.values()):
                    classification_features = [[r[f] for f in features] for r in classification_records]
                    predictions = classifier.predict_proba(classification_features)
                    for record, prediction in zip(classification_records, predictions[:,1]):
                        record[record_label] = prediction
                        record[verdict_record_label] = "Yes" if prediction >= fpr_threshold else "No"
                        
                for bench_label, (bench_ess, bench_non_ess) in benchmarks.items():
                    bench_all_records = bench_ess + bench_non_ess
                    bench_all_features = [[r[f] for f in features] for r in bench_all_records]
                    predictions = classifier.predict_proba(bench_all_features)
                    for record, prediction in zip(bench_all_records, predictions[:,1]): 
                        record[record_label] = prediction
                        record[verdict_record_label] = "Yes" if prediction >= fpr_threshold else "No"
                    
                    roc_data = roc_curve(
                        [1] * len(bench_ess)  + [0] * len(bench_non_ess),
                        predictions[:,1]
                    )
                    
                    write_auc_curve(roc_data, os.path.join(data_folder, "bench_AUC.%s.%s.png" % (bench_label, record_label)))
        
        with pd.ExcelWriter(os.path.join(data_folder, "tables.xlsx")) as excel_writer:
            # In Excel, the order of the sheets is:
            # 1) Prediction table - combined.
            # 2) Prediction table - pre.
            # 3) Predicion table - post.
            # 4) Training table
            # 5) FPR tables - classifiers X feature groups
            combined_col_config = explode_col_config(class_col_config)
            combined_col_config[-2:-2] = \
                [{"field_name": "%s-%s-verdict" % key, "csv_name": "%s - %s - final ess. verdict" % key,}
                 for key in product(cls_factory.keys(), cls_feature_groups.keys())]
            
            # TODO: can't be generalized at all!
            for ds_label, datasets in class_datasets.items():
                if len(datasets) == 2:
                    combined_records = combine_pre_post(datasets[0], datasets[1])
                    for (cls_name, grp_name) in product(cls_factory.keys(), cls_feature_groups.keys()):
                        pre_verdict_key = "pre_%s-%s-verdict" % (cls_name, grp_name)
                        post_verdict_key = "post_%s-%s-verdict" % (cls_name, grp_name)
                        combined_verdict_key = "%s-%s-verdict" % (cls_name, grp_name)
                        for r in combined_records:
                            combined_verdict = {
                                ("Yes", "Yes"): "Essential",
                                ("No", "No"): "Not essential",
                                ("Yes", "No"): "Weirdo",
                                ("No", "Yes"): "Sick"
                            }[(r[pre_verdict_key], r[post_verdict_key])]
                            r[combined_verdict_key] = combined_verdict
                    combined_df = SummaryTable.write_data_to_data_frame(combined_records, combined_col_config)
                    combined_df.to_excel(excel_writer, sheet_name="%s - combined class." % ds_label, index=False)
                else:
                    assert len(datasets) == 1
                    pred_df = SummaryTable.write_data_to_data_frame(datasets[0], class_col_config)
                    pred_df.to_excel(excel_writer, sheet_name="%s - class" % ds_label, index=False)
            
            train_df = SummaryTable.write_data_to_data_frame(train_ess_records + train_non_ess_records, train_col_config)
            train_df.to_excel(excel_writer, sheet_name="Training", index=False)
                
            for cls_name, grp_name in product(cls_factory.keys(), cls_feature_groups.keys()):
                label = "train-%s-%s" % (cls_name, grp_name)
                fprs, tprs, thresholds = roc_curve(
                    train_annotations,
                    [r[label] for r in train_all_records]
                )
                was_selected = [""] * len(fprs)
                was_selected[threshold_ixs[(cls_name, grp_name)]] = "Selected"
                fpr_df = pd.DataFrame(
                    {"FPR": fprs, "TPR": tprs, "Threshold": thresholds, "Was selected": was_selected},
                    columns=["FPR", "TPR", "Threshold", "Was selected"]
                )
                fpr_df.to_excel(excel_writer, sheet_name="FPRs - %s - %s" % (cls_name, grp_name), index=False)
                
            # Print benchmark tables
            for bench_label, (bench_ess, bench_non_ess) in benchmarks.items():
                bench_df = SummaryTable.write_data_to_data_frame(bench_ess + bench_non_ess, class_col_config)
                bench_df.to_excel(excel_writer, sheet_name="Benchmark %s - class." % bench_label, index=False)
                
                bench_all_records = bench_ess + bench_non_ess
                bench_annotations = [1] * len(bench_ess) + [0] * len(bench_non_ess)
                
                for cls_name, grp_name in product(cls_factory.keys(), cls_feature_groups.keys()):
                    label = "%s-%s" % (cls_name, grp_name)
                    fprs, tprs, thresholds = roc_curve(
                        bench_annotations,
                        [r[label] for r in bench_all_records]
                    )
                    fpr_df = pd.DataFrame(
                        {"FPR": fprs, "TPR": tprs, "Threshold": thresholds},
                        columns=["FPR", "TPR", "Threshold"]
                    )
                    fpr_df.to_excel(excel_writer, sheet_name="%s - FPRs-%s-%s" % (bench_label, cls_name, grp_name), index=False)

def write_auc_curve( (fpr, tpr, _threshold), output_file, title=None ):
    classifier_auc = auc(fpr, tpr)
    plt.plot(fpr, tpr, label="AUC = %.3f" % classifier_auc)
    plt.legend(loc="lower right")
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    
    if title is not None:
        plt.title(title)
    
    plt.savefig(output_file,
                transparent=True,
                dpi=300)
    
    plt.close()

def write_ortholog_excel(orth_df, calb_fprs_df, scer_fprs_df, spom_fprs_df, output_file):
    orth_sheet_name = "Orthologs"
    ctrl_sheet_name = "Controls"
    calb_fprs_sheet_name = "Calb FPRs"
    scer_fprs_sheet_name = "Scer FPRs"
    spom_fprs_sheet_name = "Spom FPRs"
    
    with pd.ExcelWriter(output_file) as writer:
        workbook = writer.book
        for prefix in ("Ca", "Sc", "Sp"):
            insertion_index = list(orth_df.columns.values).index("%s RF G4" % prefix) + 1
            orth_df.insert(loc=insertion_index, column='%s verdict' % prefix, value=0)
        orth_df.to_excel(writer, sheet_name=orth_sheet_name, index=False)
        ctrl_sheet = workbook.add_worksheet(ctrl_sheet_name)
        calb_fprs_df.to_excel(writer, sheet_name=calb_fprs_sheet_name, index=False)
        scer_fprs_df.to_excel(writer, sheet_name=scer_fprs_sheet_name, index=False)
        spom_fprs_df.to_excel(writer, sheet_name=spom_fprs_sheet_name, index=False)
        
        ctrl_sheet_cols = ('Organism', "FPR", "TPR", "Threshold")
        ctrl_sheet.write_row(0, 0, ctrl_sheet_cols)
        ctrl_sheet.write_row(1, 0, ("Calb", '', '', 0))
        ctrl_sheet.write_row(2, 0, ("Scer", '', '', 0))
        ctrl_sheet.write_row(3, 0, ("Spom", '', '', 0))
        ctrl_sheet.add_table(0, 0, 3, 3, {'name': 'ControlTable', "header_row": True, 'autofilter': False, "columns": [{"header": col} for col in ctrl_sheet_cols]})
        
        calb_fprs_sheet = workbook.get_worksheet_by_name(calb_fprs_sheet_name)
        scer_fprs_sheet = workbook.get_worksheet_by_name(scer_fprs_sheet_name)
        spom_fprs_sheet = workbook.get_worksheet_by_name(spom_fprs_sheet_name)
        
        calb_fprs_sheet.add_table(0, 0, calb_fprs_df.shape[0], 2, {'name': "CalbFprTable", "header_row": True, 'autofilter': False, "columns": [{"header": col} for col in list(calb_fprs_df.columns.values)]})
        scer_fprs_sheet.add_table(0, 0, scer_fprs_df.shape[0], 2, {'name': "ScerFprTable", "header_row": True, 'autofilter': False, "columns": [{"header": col} for col in list(scer_fprs_df.columns.values)]})
        spom_fprs_sheet.add_table(0, 0, spom_fprs_df.shape[0], 2, {'name': "SpomFprTable", "header_row": True, 'autofilter': False, "columns": [{"header": col} for col in list(spom_fprs_df.columns.values)]})
        
        ctrl_sheet.write_formula(1, 1, '=INDEX(CalbFprTable[FPR], MATCH(ControlTable[[#This Row],[Threshold]], CalbFprTable[Threshold], 0))')
        ctrl_sheet.write_formula(1, 2, '=INDEX(CalbFprTable[TPR], MATCH(ControlTable[[#This Row],[Threshold]], CalbFprTable[Threshold], 0))')
        ctrl_sheet.write_formula(2, 1, '=INDEX(ScerFprTable[FPR], MATCH(ControlTable[[#This Row],[Threshold]], ScerFprTable[Threshold], 0))')
        ctrl_sheet.write_formula(2, 2, '=INDEX(ScerFprTable[TPR], MATCH(ControlTable[[#This Row],[Threshold]], ScerFprTable[Threshold], 0))')
        ctrl_sheet.write_formula(3, 1, '=INDEX(SpomFprTable[FPR], MATCH(ControlTable[[#This Row],[Threshold]], SpomFprTable[Threshold], 0))')
        ctrl_sheet.write_formula(3, 2, '=INDEX(SpomFprTable[TPR], MATCH(ControlTable[[#This Row],[Threshold]], SpomFprTable[Threshold], 0))')
        
        orth_sheet = workbook.get_worksheet_by_name(orth_sheet_name)
        orth_sheet.add_table(0, 0, orth_df.shape[0], orth_df.shape[1]-1, {'name': "OrthologsTable", "header_row": True, "columns": [{"header": col} for col in list(orth_df.columns.values)]})
        for org_ix, prefix in enumerate(("Ca", "Sc", "Sp")):
            verdict_col_ix = list(orth_df.columns.values).index("%s verdict" % prefix)
            for row_ix in xrange(1, orth_df.shape[0]+1):
                orth_sheet.write_formula(row_ix, verdict_col_ix, '=IF(OrthologsTable[[#This Row],[%s RF G4]]>=\'%s\'!$D$%d, "Yes", "No")' % (prefix, ctrl_sheet_name, org_ix+2))

    
def main():
    # List of folders used for source data:
    # TODO: this should be generalized eventually and refactored with argparse,
    # but currently too much stuff is hardcoded, and in any case no major
    # changes are to be done for the paper at this point. 
    output_folder = os.path.join(Shared.get_script_dir(), "output", "predictions")
    alb_hit_file_folder = Shared.get_dependency("albicans", "experiment data", "post evo", "q20m2")
    all_cer_track_files = glob.glob(Shared.get_dependency("Kornmann", "*WildType*.wig"))
    pombe_hit_file = Shared.get_dependency("pombe", "hermes_hits.csv")
    
    # TODO: should not be hard-coded. Consider using DomainFigures to create the figures
    # as needed, instead of copying them from somewhere else.
    calb_source_figure_folder = "/Users/bermanlab/OneDrive2/OneDrive/Tn Paper/All gene figures/Calb"
    calb_source_figure_list = os.listdir(calb_source_figure_folder)
    
    scer_source_figure_folder = "/Users/bermanlab/OneDrive2/OneDrive/Tn Paper/All gene figures/Scer"
    scer_source_figure_list = os.listdir(scer_source_figure_folder)
    
    spom_source_figure_folder = "/Users/bermanlab/OneDrive2/OneDrive/Tn Paper/All gene figures/Spom"
    spom_source_figure_list = os.listdir(spom_source_figure_folder)
    
    # Set up needed infrastructure
    # The feature DBs:
    alb_db = Organisms.alb.feature_db
    cer_db = Organisms.cer.feature_db
    pom_db = Organisms.pom.feature_db
    
    # Read cerevisiae hit data:
    all_cer_tracks = [SummaryTable.get_hits_from_wig(fname) for fname in all_cer_track_files]
    cer_wt_combined = sum(all_cer_tracks, [])
    
    cer_wt_combined_analyzed = SummaryTable.analyze_hits(cer_wt_combined, cer_db, 10000).values()
    cer_wt_combined_benchmark_records = list( only_orfs(no_dubious_features(no_depleted_hit_features(cer_wt_combined_analyzed))) )
    ignored_cer_genes = Organisms.cer.ignored_features
    cer_wt_combined_benchmark_records = [r for r in cer_wt_combined_benchmark_records if r["feature"].standard_name not in ignored_cer_genes]
    
    for record in cer_wt_combined_analyzed:
        rec_name = record["feature"].standard_name
        record["ground_truth"] = "Yes" if rec_name in Organisms.cer.literature_essentials else "No" if rec_name in Organisms.cer.literature_non_essentials else ""
        record["has_paralog"] = "Yes" if rec_name in Organisms.cer.genes_with_paralogs else ""

    # Read albicans hit data:
    rdf = 1 # Read depth filter used
    post_hits = SummaryTable.read_hit_files(glob.glob(os.path.join(alb_hit_file_folder, '*03*')) +
                                            glob.glob(os.path.join(alb_hit_file_folder, '*07*')) +
                                            glob.glob(os.path.join(alb_hit_file_folder, '*11*')),
                                            rdf)
    post_analyzed = SummaryTable.analyze_hits(sum(post_hits, []), alb_db)
    
    ignored_alb_genes = Organisms.alb.ignored_features
    ignored_alb_genes |= set(f.standard_name for f in alb_db.get_all_features() if not f.is_orf or "dubious" in f.type.lower())
    for key in post_analyzed.keys():
        if key in ignored_alb_genes:
            del post_analyzed[key]
    
    SummaryTable.enrich_alb_records(post_analyzed.values())
    
    # These FPs and FNs were manually curated by us:
    alb_lit_fps = set(["C5_02260C_A", "C4_07200C_A", "C7_01840W_A", "C3_07480W_A", "C1_09370W_A", "C1_06900C_A", "C4_04930C_A", "C1_14190C_A", "C6_00320C_A", "C5_05310W_A", "C1_11280W_A", "C2_05140W_A", "C4_05180C_A", "C1_14470W_A"])
    alb_lit_fns = set(["C4_00610W_A", "C2_10210C_A", "C1_03600W_A", "C7_00700W_A", "C4_04730W_A", "CR_01740W_A", "CR_00710C_A", "C2_04760W_A", "C1_04380W_A", "C4_03440C_A", "C1_04090C_A", "C4_04180C_A", "C1_10860C_A", "C5_02900W_A", "CR_06420W_A", "C2_09370C_A", "CR_03240C_A", "C1_04330W_A", "C7_00890C_A", "C2_07100W_A", "C2_04220C_A", "C6_02840C_A", "C5_01720C_A", "C1_02230W_A", "C1_09870W_A", "C2_06540C_A", "C1_01790W_A", "C3_07550C_A", "C5_04600C_A", "C3_02960C_A", "C5_05190W_A", "CR_03430W_A", "C1_12510W_A", "CR_10140W_A", "CR_05620C_A", "C1_00700W_A", "C4_01190W_A", "C1_01490W_A", "C4_00130W_A", "C7_04230W_A", "C1_10210C_A", "CR_05030W_A", "C1_11400C_A", "C1_00060W_A", "CR_07580C_A", "C7_03940C_A", "CR_06640C_A", "C4_04090C_A", "C4_04850C_A", "C1_06230C_A", "C1_07970C_A", "C2_03180C_A"])
    
    alb_orths_essential_in_sp_and_sc = SummaryTable.get_calb_ess_in_sc()[0] & SummaryTable.get_calb_ess_in_sp()[0]
    alb_deleted_genes = SummaryTable.get_homann_deletions() | SummaryTable.get_noble_deletions() | SummaryTable.get_sanglard_deletions()
    
    # These are "systematic", but not manually filtered.
    alb_lit_ess = alb_orths_essential_in_sp_and_sc - alb_deleted_genes
    alb_lit_non_ess = alb_deleted_genes - alb_orths_essential_in_sp_and_sc
    
    # Used for training:
    filtered_alb_lit_ess = alb_lit_ess - alb_lit_fps
    filtered_alb_lit_non_ess = alb_lit_non_ess - alb_lit_fns
    
    # TODO: figure out the Roemer thing.
    roemer_ess, roemer_non_ess, omeara_ess, omeara_non_ess = \
        SummaryTable.get_grace_essentials()
    
    # TODO: generalize for cerevisiae and pombe?
    for record in post_analyzed.values():
        rec_name = record["feature"].standard_name
        record["filtered_in_training"] = "FP" if rec_name in alb_lit_fps else "FN" if rec_name in alb_lit_fns else ""
        record["ground_truth"] = "Yes" if rec_name in alb_lit_ess else "No" if rec_name in alb_lit_non_ess else ""
    
    # Read the pombe hit data:
    pom_hits = SummaryTable.read_pombe_hit_file(pombe_hit_file, rdf)
    ignored_pom_genes = Organisms.pom.ignored_features
    pom_records = SummaryTable.analyze_hits(pom_hits, pom_db)
    pom_records = {n: f for n, f in pom_records.items() if n not in ignored_pom_genes}
    pom_ess_names = Organisms.pom.literature_essentials
    pom_non_ess_names = Organisms.pom.literature_non_essentials
    pom_ess_names -= ignored_pom_genes
    pom_non_ess_names -= ignored_pom_genes
    pom_ess_records = [pom_records[n] for n in pom_ess_names if n in pom_records]
    pom_non_ess_records = [pom_records[n] for n in pom_non_ess_names if n in pom_records]
    for r in pom_records.values():
        std_name = r["feature"].standard_name 
        if std_name in pom_ess_names:
            r["ground_truth"] = "Yes"
        elif std_name in pom_non_ess_names:
            r["ground_truth"] = "No"
        else:
            r["ground_truth"] = ""
        
        r["has_paralog"] = "Yes" if std_name in Organisms.pom.genes_with_paralogs else ""
        
    # Set up the classifer:
    
    cols_config = COLS_CONFIG
    
    # Define the classification parameters:
    feature_groups = OrderedDict((
#         ("G3", ("neighborhood_index", "length", "hits", "freedom_index", "upstream_hits_100")),
        ("G4", ("neighborhood_index", "length", "hits", "reads", "freedom_index", "upstream_hits_100")),
#         ("G5", ("neighborhood_index", "length", "hits", "reads_ni", "freedom_index", "upstream_hits_100")),
    ))
    
    classifier_factory = {
#         "LR": lambda: LogisticRegression(),
        "RF": lambda: RandomForestClassifier(n_estimators=100, random_state=0)
    }
            
    pom_db = GenomicFeatures.default_pom_db()
    spom_core_ess = set(f.standard_name for f in
                        (pom_db.get_feature_by_name(r["pombe_ortholog"]) for r in post_analyzed.values()
                         if r["essential_in_pombe"] == r["essential_in_cerevisiae"] == "Yes")
                        if f) - Organisms.pom.ignored_features
    spom_core_training_ess = map(dict, (r for n, r in pom_records.items() if n in spom_core_ess))
    spom_core_non_ess = set(f.standard_name for f in
                            (pom_db.get_feature_by_name(r["pombe_ortholog"]) for r in post_analyzed.values()
                             if r["essential_in_pombe"] == r["essential_in_cerevisiae"] == "No")
                            if f) - Organisms.pom.ignored_features
    spom_core_training_non_ess = map(dict, (r for n, r in pom_records.items() if n in spom_core_non_ess))
    scer_core_ess = set(f.standard_name for f in
                        (cer_db.get_feature_by_name(list(r["feature"].cerevisiae_orthologs)[0]) for r in post_analyzed.values()
                         if r["essential_in_pombe"] == r["essential_in_cerevisiae"] == "Yes")
                        if f) - Organisms.cer.ignored_features
    scer_core_training_ess = map(dict, (r for r in cer_wt_combined_analyzed if r["feature"].standard_name in scer_core_ess))
    scer_core_non_ess = set(f.standard_name for f in
                            (cer_db.get_feature_by_name(list(r["feature"].cerevisiae_orthologs)[0]) for r in post_analyzed.values()
                             if r["essential_in_pombe"] == r["essential_in_cerevisiae"] == "No")
                            if f) - Organisms.cer.ignored_features
    scer_core_training_non_ess = map(dict, (r for r in cer_wt_combined_analyzed if r["feature"].standard_name in scer_core_non_ess))
    calb_core_training_ess = map(dict, (r for r in post_analyzed.itervalues()
                                        if r["essential_in_cerevisiae"] == "Yes" and
                                        r["essential_in_pombe"] == "Yes"))
    calb_core_training_non_ess = map(dict, (r for r in post_analyzed.itervalues()
                                            if r["essential_in_cerevisiae"] == "No" and
                                            r["essential_in_pombe"] == "No"))
    
    
    copy_calb_class_dataset = lambda: map(dict, post_analyzed.itervalues())
    copy_scer_com_class_dataset = lambda: map(dict, cer_wt_combined_benchmark_records)
    copy_spom_class_dataset = lambda: map(dict, (r for k, r in pom_records.iteritems() if k not in ignored_pom_genes))
    
    calb_ortholog_class_dataset = copy_calb_class_dataset()
    scer_ortholog_class_dataset = copy_scer_com_class_dataset()
    spom_ortholog_class_dataset = copy_spom_class_dataset()
    
    copy_omeara_bench_ess = lambda: [dict(r) for r in post_analyzed.values() if r["feature"].standard_name in omeara_ess]
    copy_omeara_bench_non_ess = lambda: [dict(r) for r in post_analyzed.values() if r["feature"].standard_name in omeara_non_ess]
    
    copy_calb_bench_ess = lambda:  map(dict, (r for n, r in post_analyzed.iteritems() if n in filtered_alb_lit_ess))
    copy_calb_bench_non_ess = lambda: map(dict, (r for n, r in post_analyzed.iteritems() if n in filtered_alb_lit_non_ess))
    
    training_calb_manually_curated_ess = [dict(post_analyzed[n]) for n in filtered_alb_lit_ess if n in post_analyzed]
    training_calb_manually_curated_non_ess = [dict(post_analyzed[n]) for n in filtered_alb_lit_non_ess if n in post_analyzed]
    
    # These FPs and FNs were manually curated by us:
    scer_training_fps = set(["S000005081", "S000000510", "S000003278", "S000005250", "S000002448", "S000006163"])
    scer_training_fns = set(["S000000370", "S000005217", "S000003026", "S000005194", "S000000175", "S000002583", "S000003092", "S000001390", "S000003850", "S000000223", "S000001915", "S000004391", "S000001295", "S000003716", "S000003148", "S000003073", "S000002856", "S000003063", "S000003440", "S000005710", "S000000393", "S000002248", "S000005850", "S000000363", "S000006169", "S000004388", "S000006335", "S000002653", "S000005024", "S000001537", "S000003293", "S000000070"])
    
    for r in cer_wt_combined_benchmark_records:
        name = r["feature"].standard_name
        r["filtered_in_training"] = "FP" if name in scer_training_fps else "FN" if name in scer_training_fns else ""
    
    # These FPs and FNs were manually curated by us:
    spom_training_fps = set(["SPCC330.10", "SPBC14C8.14c", "SPBC725.17c", "SPBC25H2.04c", "SPAC144.07c", "SPAC31A2.05c", "SPBC146.05c", "SPAC144.18", "SPCC18.12c", "SPCC1450.10c", "SPBC25H2.06c", "SPAC22E12.10c", "SPAC4D7.12c", "SPBC30D10.02", "SPBC3B9.12", "SPBC16G5.10", "SPAC3G9.12", "SPCC16C4.08c", "SPAC17A5.13", "SPAC1F5.02", "SPBC685.05", "SPBC21.06c", "SPAC806.02c", "SPAC16E8.15", "SPBC19F8.07", "SPCC63.10c", "SPBC211.01", "SPBC2G2.04c", "SPAC31G5.08", "SPAC1006.02", "SPBC14F5.04c", "SPBC9B6.10", "SPBC3B9.21", "SPAC16A10.06c", "SPBC21C3.10c", "SPAC56E4.02c", "SPAC57A7.11", "SPAC144.08", "SPBC3B8.01c", "SPAC3G9.16c", "SPBC16D10.10", "SPBC36.12c", "SPAC12G12.04", "SPAC6F12.05c", "SPAC222.03c", "SPBC1734.03"])
    spom_training_fns = set(["SPAC13A11.03", "SPAC1F7.07c", "SPCC1442.03", "SPBC646.13", "SPAC19G12.08", "SPAC1F12.07", "SPCC191.02c", "SPBC405.04c", "SPAC890.06", "SPAC1B2.03c", "SPBC119.05c", "SPAP27G11.05c", "SPBC23G7.08c", "SPBC19F8.03c", "SPAC977.17", "SPAC30D11.10", "SPBC530.13", "SPAC18G6.04c", "SPBC56F2.10c", "SPCC736.06", "SPCC622.18", "SPCC737.02c", "SPBC29A3.01", "SPAC6B12.15", "SPBC2G2.01c", "SPBC119.06", "SPBC4B4.03", "SPAC328.02", "SPAC664.02c", "SPAC20H4.07", "SPAC8C9.03", "SPAC227.01c", "SPAC25H1.07", "SPAC16E8.13", "SPCC594.05c", "SPBC776.03", "SPAC1952.09c", "SPAC2F7.07c", "SPAC1F3.10c", "SPAC22F3.09c", "SPCC830.06", "SPCC594.06c", "SPAC4G8.10", "SPAC4F8.01", "SPCC16C4.09", "SPAC20G8.10c", "SPBC146.13c", "SPAC1B3.07c", "SPBC3F6.05", "SPBC215.05", "SPAC17G6.04c", "SPBC1706.03", "SPBP16F5.07", "SPAC644.14c", "SPCC61.02", "SPAC23C11.08", "SPBC1778.06c", "SPCC338.14", "SPCC18B5.03", "SPAC15A10.03c", "SPBC32F12.01c", "SPBC887.10"])
    
    for n, r in pom_records.iteritems():
        r["filtered_in_training"] = "FP" if n in spom_training_fps else "FN" if n in spom_training_fns else ""
    
    filtered_scer_core_training_ess = [dict(r) for r in scer_core_training_ess if r["feature"].standard_name not in scer_training_fps]
    filtered_scer_core_training_non_ess = [dict(r) for r in scer_core_training_non_ess if r["feature"].standard_name not in scer_training_fns]
    
    filtered_spom_core_training_ess = [dict(r) for r in spom_core_training_ess if r["feature"].standard_name not in spom_training_fps]
    filtered_spom_core_training_non_ess = [dict(r) for r in spom_core_training_non_ess if r["feature"].standard_name not in spom_training_fns]
    
    copy_scer_com_bench_ess = lambda: [dict(r) for r in cer_wt_combined_benchmark_records if r["feature"].standard_name in (Organisms.cer.literature_essentials - scer_training_fps)]
    copy_scer_com_bench_non_ess = lambda: [dict(r) for r in cer_wt_combined_benchmark_records if r["feature"].standard_name in (Organisms.cer.literature_non_essentials - scer_training_fns)]
    
    copy_spom_bench_ess = lambda: map(dict, (r for r in pom_ess_records if r["feature"].standard_name not in spom_training_fps))
    copy_spom_bench_non_ess = lambda: map(dict, (r for r in pom_non_ess_records if r["feature"].standard_name not in spom_training_fns))
    
    data = OrderedDict((
        ("Calb training - manually curated", (
            training_calb_manually_curated_ess,
            training_calb_manually_curated_non_ess,
            OrderedDict((
                ("Calb-all", (copy_calb_class_dataset(),)),
                ("Calb-ortholog-only", (calb_ortholog_class_dataset,)),
            )),
            OrderedDict((
                ("Scer-com", (
                    copy_scer_com_bench_ess(),
                    copy_scer_com_bench_non_ess(),
                )),
                ("Spom", (
                    copy_spom_bench_ess(),
                    copy_spom_bench_non_ess(),
                )),
                ("Calb", (
                    copy_calb_bench_ess(),
                    copy_calb_bench_non_ess()
                )),
            ))
        )),
        ("Scer training - filtered", (
            filtered_scer_core_training_ess,
            filtered_scer_core_training_non_ess,
            OrderedDict((
                ("Calb", (copy_calb_class_dataset(),)),
                ("Scer-com", (scer_ortholog_class_dataset,)),
                ("Spom", (copy_spom_class_dataset(),)),
            )),
            OrderedDict((
                ("Scer-com", (
                    copy_scer_com_bench_ess(),
                    copy_scer_com_bench_non_ess(),
                )),
                ("Spom", (
                    copy_spom_bench_ess(),
                    copy_spom_bench_non_ess(),
                )),
                ("Calb", (
                    copy_calb_bench_ess(),
                    copy_calb_bench_non_ess()
                )),
            ))
        )),
        ("Spom training - filtered", (
            filtered_spom_core_training_ess,
            filtered_spom_core_training_non_ess,
            OrderedDict((
                ("Calb", (copy_calb_class_dataset(),)),
                ("Scer-com", (copy_scer_com_class_dataset(),)),
                ("Spom", (spom_ortholog_class_dataset,)),
            )),
            OrderedDict((
                ("Scer-com", (
                    copy_scer_com_bench_ess(),
                    copy_scer_com_bench_non_ess(),
                )),
                ("Spom", (
                    copy_spom_bench_ess(),
                    copy_spom_bench_non_ess(),
                )),
                ("Calb", (
                    copy_calb_bench_ess(),
                    copy_calb_bench_non_ess()
                )),
            ))
        )),
    ))

    run_pipeline(classifier_factory, feature_groups, data, cols_config, cols_config, output_folder, 0.1)
    
    # TODO: what does this do? Add a comment.

    inter_col_config = list(cols_config)
    inter_col_config[-2:-2] = [
        {"field_name": "ground_truth", "csv_name": "Ess. ground truth"},
        {"field_name": "RF-G4", "csv_name": "RF - G4" , "format": "%.3f"},
        {"field_name": "RF-G4-verdict", "csv_name": "RF - G4 - ess. verdict"}
    ]
    
    # Add intersection folders using tn screen
    # Add literature vs. predictor in all three (use O'Meara as literature,
    # also "List of possibly ess genes from Aaron Mitchell.xlsx" in Tn paper.
    
    score_to_use = "RF-G4"
    verdict_to_use = "%s-verdict" % score_to_use
    train_score_to_use = "train-%s" % score_to_use
    
    # Draw Venns:
    calb_records = {r["feature"].standard_name: r for r in calb_ortholog_class_dataset}.items()
    calb_records_only = [p[1] for p in calb_records]
    orthologs = set(k for k, r in post_analyzed.items() if r["pombe_ortholog"] and r["feature"].cerevisiae_orthologs)
    draw_venn(
        "Essential orthologs",
        os.path.join(output_folder, "orthologs_calb_ess.png"),
        [set(k for k, r in calb_records if r[verdict_to_use] == "Yes") & orthologs,
         set(k for k, r in calb_records if r["essential_in_cerevisiae"] == "Yes") & orthologs,
         set(k for k, r in calb_records if r["essential_in_pombe"] == "Yes") & orthologs],
        ["Calb", "Scer", "Spom"]
    )
    
    draw_venn(
        "Non-Essential orthologs",
        os.path.join(output_folder, "orthologs_calb_non_ess.png"),
        [set(k for k, r in calb_records if r[verdict_to_use] == "No") & orthologs,
         set(k for k, r in calb_records if r["essential_in_cerevisiae"] == "No") & orthologs,
         set(k for k, r in calb_records if r["essential_in_pombe"] == "No") & orthologs],
        ["Calb", "Scer", "Spom"]
    )
    
    calb_essentials = set(k for k, r in post_analyzed.items() if r["essential_in_albicans_grace_omeara"] and r["essential_in_mitchell"])
    draw_venn(
        "Essential in Calb",
        os.path.join(output_folder, "tn_omeara_mitchell_calb_ess.png"),
        [set(k for k, r in calb_records if r[verdict_to_use] == "Yes") & calb_essentials,
         set(k for k, r in calb_records if r["essential_in_albicans_grace_omeara"] == "Yes") & calb_essentials,
         set(k for k, r in calb_records if r["essential_in_mitchell"] == "Yes") & calb_essentials,],
        ["Tn", "O'Meara", "Mitchell"]
    )
    
    draw_venn(
        "Cmp. with Mitchell ess.",
        os.path.join(output_folder, "tn_mitchell_calb_ess.png"),
        [set(k for k, r in calb_records if r[verdict_to_use] == "Yes" and r["essential_in_mitchell"]),
         set(k for k, r in calb_records if r["essential_in_mitchell"] == "Yes")],
        ["Tn", "Mitchell"]
    )
    
    draw_venn(
        "Cmp. with O'Meara ess.",
        os.path.join(output_folder, "tn_omeara_calb_ess.png"),
        [set(k for k, r in calb_records if r[verdict_to_use] == "Yes" and r["essential_in_albicans_grace_omeara"]),
         set(k for k, r in calb_records if r["essential_in_albicans_grace_omeara"] == "Yes")],
        ["Tn", "O'Meara"]
    )
    
    draw_venn(
        "Cmp. with O'Meara non-ess.",
        os.path.join(output_folder, "tn_omeara_calb_non_ess.png"),
        [set(k for k, r in calb_records if r[verdict_to_use] == "No" and r["essential_in_albicans_grace_omeara"]),
         set(k for k, r in calb_records if r["essential_in_albicans_grace_omeara"] == "No")],
        ["Tn", "O'Meara"]
    )
    
    def dump_records(
            records,
            source_figure_folder,
            source_figure_list,
            inter_out_folder,
            add_score=False
        ):
        target_figure_folder = os.path.join(inter_out_folder, "figures")
        if not os.path.exists(target_figure_folder):
            Shared.make_dir(target_figure_folder)
        
        for record in records:
            feature = record["feature"]
            gene_name = getattr(feature, "feature_name", feature.standard_name)
            for source_figure in source_figure_list:
                if gene_name in source_figure:
                    break
            else:
                gene_name = None
            
            # TODO: the score is hard-coded :(
            dest = target_figure_folder if not add_score else \
                os.path.join(target_figure_folder, "%.3f_%s" % (record["train-"+score_to_use], source_figure))
            shutil.copy(os.path.join(source_figure_folder, source_figure), dest)
        
        with pd.ExcelWriter(os.path.join(inter_out_folder, "genes.xlsx")) as excel_writer:
            sheet = SummaryTable.write_data_to_data_frame(records, inter_col_config)
            sheet.to_excel(excel_writer, sheet_name="Main", index=False)
    
    for ca_status, sc_status, sp_status in product(("Yes", "No"), repeat=3):
        inter_out_folder = os.path.join(output_folder, "Ca-%s Sc-%s Sp-%s" % (ca_status, sc_status, sp_status))
        if not os.path.exists(inter_out_folder):
            os.mkdir(inter_out_folder)
        
        records = [r for r in calb_records_only if
                   r["essential_in_cerevisiae"] == sc_status and
                   r["essential_in_pombe"] == sp_status and
                   r[verdict_to_use] == ca_status]
        
        dump_records(records, calb_source_figure_folder, calb_source_figure_list, inter_out_folder)
        
    for tn_status, omeara_status in product(("Yes", "No"), repeat=2):
        inter_out_folder = os.path.join(output_folder, "Mitchell-Yes Calb-Tn-%s OMeara-%s" % (tn_status, omeara_status))
        if not os.path.exists(inter_out_folder):
            os.mkdir(inter_out_folder)
        
        records = [r for r in calb_records_only if
                   r["essential_in_mitchell"] == "Yes" and
                   r["essential_in_albicans_grace_omeara"] == omeara_status and
                   r[verdict_to_use] == tn_status]
        
        dump_records(records, calb_source_figure_folder, calb_source_figure_list, inter_out_folder)
    
    inter_col_config[-2:-2] = [
        {"field_name": "train-RF-G4", "csv_name": "RF - G4" , "format": "%.3f"},
        {"field_name": "train-RF-G4-verdict", "csv_name": "RF - G4 - ess. verdict"}
    ]
    
    # Print out the conflicting annotations in the benchmark datasets:
    for training_label, dataset in data.iteritems():
        training_ess, training_non_ess = dataset[:2]        
        
        dump_records(
            training_ess,
            {"Calb": calb_source_figure_folder, "Scer": scer_source_figure_folder, "Spom": spom_source_figure_folder}[training_label[:4]],
            {"Calb": calb_source_figure_list, "Scer": scer_source_figure_list, "Spom": spom_source_figure_list}[training_label[:4]],
            os.path.join(output_folder, "classification - %s" % training_label, "Training essentials"),
            add_score=True
        )
        dump_records(
            training_non_ess,
            {"Calb": calb_source_figure_folder, "Scer": scer_source_figure_folder, "Spom": spom_source_figure_folder}[training_label[:4]],
            {"Calb": calb_source_figure_list, "Scer": scer_source_figure_list, "Spom": spom_source_figure_list}[training_label[:4]],
            os.path.join(output_folder, "classification - %s" % training_label, "Training non-essentials"),
            add_score=True
        )
        
    # Ortholog megatable:
    
    # TODO: we already had this indexed... bah...
    index_dataset = lambda dataset: {r["feature"].standard_name: r for r in dataset}
    ortholog_records = []
    indexed_calb = index_dataset(calb_ortholog_class_dataset)
    indexed_scer = index_dataset(scer_ortholog_class_dataset)
    indexed_spom = index_dataset(spom_ortholog_class_dataset)
    
    # In this code we only generated the orthologs. Currently commented out
    # because what we want is a left-join of the Calb records with Scer and Spom.
#     orthologs = Organisms.get_all_orthologs()
#     for (calb_orth, scer_orth, spom_orth) in orthologs:
#         calb_record = indexed_calb.get(calb_orth.standard_name)
#         scer_record = indexed_scer.get(scer_orth.standard_name)
#         spom_record = indexed_spom.get(spom_orth.standard_name)
#         
#         if not calb_record or not scer_record or not spom_record:
#             continue
#         
#         combined_record = {}
#         for prefix, record in zip(("ca", "sc", "sp"), (calb_record, scer_record, spom_record)):
#             for key, value in record.iteritems():
#                 combined_record["%s_%s" % (prefix, key)] = value
#         
#         ortholog_records.append(combined_record)

    # The left-join implementation
    orthologs = [Organisms.get_orths_by_name(n) for n in indexed_calb.keys()]
    sc_rec_items_to_add = zip(indexed_scer.values()[0].keys(), repeat(""))
    sp_rec_items_to_add = zip(indexed_spom.values()[0].keys(), repeat(""))
    for (calb_orth, scer_orth, spom_orth) in orthologs:
        calb_record = indexed_calb.get(calb_orth.standard_name)
        scer_record = indexed_scer.get(scer_orth.standard_name if scer_orth else "")
        spom_record = indexed_spom.get(spom_orth.standard_name if spom_orth else "")
        
        combined_record = {}
        for prefix, record in zip(("ca", "sc", "sp"), (calb_record, scer_record, spom_record)):
            if record is not None:
                items_to_add = record.iteritems()
            else:
                if prefix == "sc":
                    items_to_add = sc_rec_items_to_add
                elif prefix == "sp":
                    items_to_add = sp_rec_items_to_add
            
            for key, value in items_to_add:
                combined_record["%s_%s" % (prefix, key)] = value
         
        ortholog_records.append(combined_record)
    
    # We do want to left-join Calb and the other 
    
    orth_cols_config = ORTH_COLS_CONFIG
    
    def get_roc_curve_from_data(data, name):
        training_ess = data[name][0]
        training_non_ess = data[name][1]
        return roc_curve(
            [1]*len(training_ess) + [0]*len(training_non_ess),
            [r[train_score_to_use] for r in training_ess] + [r[train_score_to_use] for r in training_non_ess]  
        )
    
    write_ortholog_excel(
        SummaryTable.write_data_to_data_frame(ortholog_records, orth_cols_config),
        pd.DataFrame(
            dict(zip(("FPR", "TPR", "Threshold"), get_roc_curve_from_data(data, "Calb training - manually curated"))),
            columns=("FPR", "TPR", "Threshold")
        ),
        pd.DataFrame(
            dict(zip(("FPR", "TPR", "Threshold"), get_roc_curve_from_data(data, "Scer training - filtered"))),
            columns=("FPR", "TPR", "Threshold")
        ),
        pd.DataFrame(
            dict(zip(("FPR", "TPR", "Threshold"), get_roc_curve_from_data(data, "Spom training - filtered"))),
            columns=("FPR", "TPR", "Threshold")
        ),
        os.path.join(output_folder, "orthologs_combined_ex.xlsx")
    )



COLS_CONFIG = [
        {
            "field_name": "feature",
            "csv_name": "Standard name",
            "format": lambda f: f.standard_name
        },
        
        {
            "field_name": "feature",
            "csv_name": "Common name",
            "format": lambda f: f.common_name
        },
        
        {
            "field_name": "feature",
            "csv_name": "Sc ortholog",
            "format": lambda f: ','.join(f.cerevisiae_orthologs) if hasattr(f, 'cerevisiae_orthologs') else ""
        },
        
        {
            "field_name": "feature",
            "csv_name": "Sc std name",
            "format": lambda f: Organisms.cer.feature_db.get_feature_by_name(list(f.cerevisiae_orthologs)[0]).standard_name if len(getattr(f, 'cerevisiae_orthologs', [])) > 0 else ""
        },
        
        {
            "field_name": "cer_fitness",
            "csv_name": "Sc fitness",
        },
        
        {
            "field_name": "essential_in_cerevisiae",
            "csv_name": "Essential in Sc",
        },
        
        {
            "field_name": "cer_synthetic_lethal",
            "csv_name": "SL in Sc",
        },
        
        {
            "field_name": "pombe_ortholog",
            "csv_name": "Sp ortholog",
        },
         
        {
            "field_name": "essential_in_pombe",
            "csv_name": "Essential in Sp",
        },
        
        {
            "field_name": "essential_in_albicans_grace_omeara",
            "csv_name": "Essential in O'Meara",
        },
        
        {
            "field_name": "essential_in_mitchell",
            "csv_name": "Essential in Mitchell",
        },
        
        {
            "field_name": "deleted_in_calb",
            "csv_name": "Deleted in Calb",
        },
        
        {
            "field_name": "filtered_in_training",
            "csv_name": "Filtered in training",
        },
        
        {
            "field_name": "feature",
            "csv_name": "Length",
            "format": lambda f: len(f)
        },
        
        {
            "field_name": "hits",
            "csv_name": "Hits",
            "format": "%d",
            "local": True
        },
        
        {
            "field_name": "reads",
            "csv_name": "Reads",
            "format": "%d",
            "local": True
        },
        
        {
            "field_name": "max_free_region",
            "csv_name": "Longest free interval",
            "format": "%d",
            "local": True
        },
                   
        {
            "field_name": "freedom_index",
            "csv_name": "Freedom index",
            "format": "%.2f",
            "local": True
        },
        
        {
            "field_name": "neighborhood_index",
            "csv_name": "Neighborhood index",
            "format": "%.3f",
            "local": True
        },

        {
            "field_name": "upstream_hits_100",
            "csv_name": "Upstream hits 100",
            "format": "%d",
            "local": True
        },
        
        {
            "field_name": "upstream_hits_50",
            "csv_name": "Upstream hits 50",
            "format": "%d",
            "local": True
        },
        
        {
            "field_name": "n_term_hits_100",
            "csv_name": "N-term hits 100",
            "format": "%d",
            "local": True
        },
        
        {
            "field_name": "feature",
            "csv_name": "Type",
            "format": lambda f: f.type
        },
        
        {
            "field_name": "feature",
            "csv_name": "Description",
            "format": lambda f: f.description
        },
    ]
    
ORTH_COLS_CONFIG = [
        {
            "field_name": "ca_feature",
            "csv_name": "Ca standard name",
            "format": lambda f: f.standard_name
        },
        
        {
            "field_name": "ca_feature",
            "csv_name": "Ca common name",
            "format": lambda f: f.common_name
        },
        
        {
            "field_name": "sc_feature",
            "csv_name": "Sc standard name",
            "format": lambda f: f.standard_name if f else ""
        },
        
        {
            "field_name": "sc_feature",
            "csv_name": "Sc common name",
            "format": lambda f: f.common_name  if f else "" or f.feature_name  if f else ""
        },
        
        {
            "field_name": "sp_feature",
            "csv_name": "Sp standard name",
            "format": lambda f: f.standard_name if f else ""
        },
        
        {
            "field_name": "sp_feature",
            "csv_name": "Sp common name",
            "format": lambda f: f.common_name if f else ""
        },
        
        {
            "field_name": "ca_essential_in_albicans_grace_omeara",
            "csv_name": "Essential in O'Meara",
        },
        
        {
            "field_name": "ca_essential_in_mitchell",
            "csv_name": "Essential in Mitchell",
        },
        
        {
            "field_name": "ca_deleted_in_calb",
            "csv_name": "Deleted in Calb",
        },
        
        {
            "field_name": "ca_has_paralog",
            "csv_name": "Paralog in Calb",
        },
        
        {
            "field_name": "ca_essential_in_cerevisiae",
            "csv_name": "Essential in Sc",
        },
        
        {
            "field_name": "ca_cer_synthetic_lethal",
            "csv_name": "SL in Sc",
        },
        
        {
            "field_name": "sc_has_paralog",
            "csv_name": "Paralog in Sc",
        },
        
        {
            "field_name": "ca_essential_in_pombe",
            "csv_name": "Essential in Sp",
        },
        
        {
            "field_name": "sp_has_paralog",
            "csv_name": "Paralog in Spom",
        },
        
        {
            "field_name": "ca_feature",
            "csv_name": "Ca length",
            "format": lambda f: len(f)
        },
        
        {
            "field_name": "sc_feature",
            "csv_name": "Sc length",
            "format": lambda f: len(f) if f else ""
        },
        
        {
            "field_name": "sp_feature",
            "csv_name": "Sp length",
            "format": lambda f: len(f) if f else ""
        },
        
        {
            "field_name": "ca_hits",
            "csv_name": "Ca hits",
            "format": "%d",
        },
        
        {
            "field_name": "sc_hits",
            "csv_name": "Sc hits",
            "format": "%d",
        },
        
        {
            "field_name": "sp_hits",
            "csv_name": "Sp hits",
            "format": "%d",
        },
        
        {
            "field_name": "ca_reads",
            "csv_name": "Ca reads",
            "format": "%d",
        },
        
        {
            "field_name": "sc_reads",
            "csv_name": "Sc reads",
            "format": "%d",
        },
        
        {
            "field_name": "sp_reads",
            "csv_name": "Sp reads",
            "format": "%d",
        },
        
        {
            "field_name": "ca_max_free_region",
            "csv_name": "Ca longest free interval",
            "format": "%d",
        },
        
        {
            "field_name": "sc_max_free_region",
            "csv_name": "Sc longest free interval",
            "format": "%d",
        },
        
        {
            "field_name": "sp_max_free_region",
            "csv_name": "Sp longest free interval",
            "format": "%d",
        },
                   
        {
            "field_name": "ca_freedom_index",
            "csv_name": "Ca freedom index",
            "format": "%.2f",
        },
        
        {
            "field_name": "sc_freedom_index",
            "csv_name": "Sc freedom index",
            "format": "%.2f",
        },
        
        {
            "field_name": "sp_freedom_index",
            "csv_name": "Sp freedom index",
            "format": "%.2f",
        },
        
        {
            "field_name": "ca_neighborhood_index",
            "csv_name": "Ca neighborhood index",
            "format": "%.3f"
        },
        
        {
            "field_name": "sc_neighborhood_index",
            "csv_name": "Sc neighborhood index",
            "format": "%.3f"
        },
        
        {
            "field_name": "sp_neighborhood_index",
            "csv_name": "Pp neighborhood index",
            "format": "%.3f"
        },
        
        {
            "field_name": "ca_upstream_hits_100",
            "csv_name": "Ca upstream hits 100",
            "format": "%d",
        },
        
        {
            "field_name": "sc_upstream_hits_100",
            "csv_name": "Sc upstream hits 100",
            "format": "%d",
        },
        
        {
            "field_name": "sp_upstream_hits_100",
            "csv_name": "Sp upstream hits 100",
            "format": "%d",
        },
        
        {
            "field_name": "ca_RF-G4",
            "csv_name": "Ca RF G4",
            "format": "%.3f",
        },
        
        {
            "field_name": "sc_RF-G4",
            "csv_name": "Sc RF G4",
            "format": "%.3f",
        },
        
        {
            "field_name": "sp_RF-G4",
            "csv_name": "Sp RF G4",
            "format": "%.3f",
        },
        
        {
            "field_name": "ca_feature",
            "csv_name": "Ca description",
            "format": lambda f: f.description
        },
        
        {
            "field_name": "sc_feature",
            "csv_name": "Sc description",
            "format": lambda f: f.description if f else ""
        },
        
        {
            "field_name": "sp_feature",
            "csv_name": "Sp description",
            "format": lambda f: f.description if f else ""
        },
    ]

if __name__ == "__main__":
    main()