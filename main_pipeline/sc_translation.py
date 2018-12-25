import pandas as pd

data = pd.read_csv('Both_target_seq_counts_test.csv')
trans = pd.read_csv('nomenclature_translation_SGD.csv')

data2 = data.merge(trans, on='id')
data2.to_csv('test.csv')

train = pd.read_csv('training_set_Sc.csv')

training = train.merge(data2)
training.to_csv('Both_ml.csv')

# training = pd.read_csv('/home/user/Desktop/transposon/main_pipeline/glabrata_training_set.csv')
# tn_features = pd.read_csv('/home/user/Desktop/transposon/main_pipeline/FreadsWT_target_seq_counts_test.csv')
# curated = pd.read_csv('curation.csv')
#
# final = training.merge(tn_features)
#
# print(final[final['id'].duplicated() == True])
#
# # final.to_csv('curation.csv')
#
# tr = set(training['id'])
# cu = set(curated['id'])
#
# removed = pd.DataFrame(pd.Series(list(tr.difference(cu)), name='id')).merge(training, tn_features, on='id')
#
# removed.to_csv('removed.csv')
