import pandas as pd

data = pd.read_csv('FreadI1_target_seq_counts_test.csv')
trans = pd.read_csv('nomenclature_translation_SGD.csv')

data2 = data.merge(trans, on='id')
data2.to_csv('test.csv')

train = pd.read_csv('training_set_Sc.csv')

training = train.merge(data2)
training.to_csv('text_ml.csv')