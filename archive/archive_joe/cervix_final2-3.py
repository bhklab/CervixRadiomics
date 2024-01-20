# import dependencies
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
# !pip install pymrmre
import pymrmre
import seaborn as sns
# !pip install lifelines scikit-learn
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
import pickle, json, time, datetime
# from google.colab import drive
# drive.mount('/content/drive')
# %matplotlib inline
# set working directory
# %cd drive/My\ Drive/notebooks/Cervix/
# %cd /Users/josephmarsilla/Google\ Drive/notebooks/Cervix/

def save_obj_json(obj, name):
    with open(name, 'w') as fp:
        json.dump(obj, fp)

def load_obj(name):
    with open(name, 'rb') as f:
        return pickle.load(f)

def save_obj(obj, name):
    with open(name, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def normalize(train_set_, test_set_, cohort_norm=False,
              ignore = ['cohort', 'figo2018', 'event', 'time', 'shape_VoxelVolume']):

  if cohort_norm:

      normalized_trainset = {}
      normalized_testset = {}

      for i, val in enumerate(a):
          b = []
          c = []
          if val not in ignore:

              for j in range(0,2):
                  # normalize by cohort...
                  me = train_set_[train_set_['cohort']==j][val].mean()
                  stds = train_set_[train_set_['cohort']==j][val].std()
                  b += list((train_set_[train_set_['cohort']==j][val] - me)/stds)
                  c += list((test_set_[test_set_['cohort']==j][val] - me)/stds)

              normalized_trainset[val] = b
              normalized_testset[val] = c

      normalized_trainset['figo2018'] = train_set_['figo2018']
      normalized_testset['figo2018'] = test_set_['figo2018']
      normalized_trainset['cohort'] = train_set_['cohort']
      normalized_testset['cohort'] = test_set_['cohort']
      train_set_ = pd.DataFrame.from_dict(normalized_trainset)
      test_set_ = pd.DataFrame.from_dict(normalized_testset)
      train_set_ = train_set_.astype('double')
      test_set_ = test_set_.astype('double')
  else:
      # ignore = [] # 'cohort', 'figo2018'
      for i, val in enumerate(a):
          if val not in ignore:
              me = train_set_[val].mean()
              stds = train_set_[val].std()
              train_set_[val] = (train_set_[val] - me)/stds
              test_set_[val]  = (test_set_[val] - me)/stds

  return train_set_, test_set_

def get_signatures(train_set_, test_set_, train_set_target_,  test_set_target_):
      lab = {}

      # lab[1] = ['shape_VoxelVolume']
      # a = list(train_mrmr)

      for l in range(2,7):
        doc = {}
        avg_c1, avg_c2 = 0.0, 0.0

        print(f'Starting {l}') # .drop(columns=['shape_VoxelVolume'])

        s1 = pymrmre.mrmr_ensemble(features=train_set_.drop(columns=['shape_VoxelVolume', 'figo2018']), targets = train_set_target_.drop(columns=['time','figo2018']),
                                    solution_length=l, solution_count=1, fixed_features=['cohort'])

        doc['features'] = list(s1.iloc[0])[0][1:]
        doc['num_features'] = l - 2 # in addition to volume...

        exclude = ['cohort', 'time', 'event','shape_VoxelVolume', 'figo2018']

        # for each value ...
        train_mrmr = train_set_[s1.iloc[0][0]]
        # what we'll be training on...
        train_mrmr['event'] = train_set_target_['event']
        train_mrmr['time'] = train_set_target_['time']

        test_mrmr = test_set_[s1.iloc[0][0]]
        test_mrmr['event'] = test_set_target_['event']
        test_mrmr['time'] = test_set_target_['time']

        b = train_mrmr.copy()
        btest = test_mrmr.copy()

        # fit cox model with chosen features
        cph1 = CoxPHFitter() #penalizer=0.01
        cph1.fit(train_mrmr, duration_col='time', event_col='event', show_progress=False) # step_size=0.25
        # doc['train_ci'] = cph1.concordance_index_
        coef = cph1.params_
        doc['cox_coef'] = coef

        a = list(train_mrmr)
        choice = ['shape_VoxelVolume']

        for c, val in enumerate(a):

          # if val in choice:
          #     b[val] *= coef[val]
          #     btest[val] *= coef[val]

          if val not in exclude:
              v =  train_mrmr[val]*coef[val]
              v_test =  test_mrmr[val]*coef[val]
              # calculates signature from all the features ...
              try:
                  b['Signature'] += v
                  btest['Signature'] += v_test
              except Exception:
                  b['Signature'] = v
                  btest['Signature'] = v_test

        # if l == 2:
        #   # choose = 'shape_VoxelVolume'
        #   choose = 'figo2018'
        # else:
        choose = "Signature"

        # volume only analysis...

        # add dicotomized based on median of training/test signature values...
        # should divide everything quite evenly in half ...
        b[choose].median()
        # b["vol_dico"] = (b['shape_VoxelVolume'] > b['shape_VoxelVolume'].median())
        b["high_signature"] = (b[choose] > b[choose].median())
        b["low_signature"] = (b[choose] < b[choose].median())
        # btest["vol_dico"] = (btest['shape_VoxelVolume'] > b['shape_VoxelVolume'].median())
        btest["high_signature"] = (btest[choose] > b[choose].median())
        btest["low_signature"] = (btest[choose] < b[choose].median())

        # run with signature for l-2 features ...
        # if l == 2:
        #   choosen = ['event','time', 'cohort', 'figo2018']
        #   # choosen = ['event','time', 'figo2018']
        #   # initially both are the same...
        #   dico_choosen = ['event', 'time', 'cohort', 'figo2018']
        # else:
        choosen = ['event','time', 'cohort', 'Signature']
        dico_choosen = ['event','time', 'cohort', "high_signature"]
        # choosen = ['event','time', 'figo2018', 'Signature']
        # dico_choosen = choosen = ['event','time', 'figo2018', "high_signature"]

        # choosen = ['event','time', 'Signature']
        # dico_choosen = choosen = ['event','time', "high_signature"]
        # add cohort here...

        cphf =  CoxPHFitter(penalizer=0.01) #penalizer=0.01
        # if l == 2:
        cphf.fit(b[choosen], duration_col='time', event_col='event', show_progress=False) # step_size=0.1
        doc['train_ci'] = cphf.concordance_index_
        print('Train CI: ', cphf.concordance_index_)
        test_signature = -cphf.predict_partial_hazard(btest[choosen]).values
        # print(test_signature)
        c1 = concordance_index(btest['time'], test_signature, btest['event'])
        doc['test_ci'] = c1
        print('Test CI: ', c1)

        # run with dicotomized signature for l - 2 features...
        # run with signature for l-2 features ...
        cphd =  CoxPHFitter(penalizer=0.01) #penalizer=0.01
        # if l == 2:
        cphd.fit(b[dico_choosen], duration_col='time', event_col='event', show_progress=False) # step_size=0.1
        doc['dicho_train_ci'] = cphd.concordance_index_
        print('Dichotomized train CI: ', cphd.concordance_index_)
        test_signature2 = -cphd.predict_partial_hazard(btest[dico_choosen]).values
        # print(test_signature)
        c2 = concordance_index(btest['time'], test_signature2, btest['event'])
        doc['dicho_test_ci'] = c2
        print('Dichotomized test CI: ',c2)

        lab[l-1] = doc
        print(f'Done {l}')

      return lab

# import data
select = pd.read_csv('/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/Features_with_all_ICC_above_0.5.csv')
combined_feat = pd.read_csv('/Users/josephmarsilla/Google Drive/notebooks/Cervix/Combinied_data_for_Joe_2020_06_16.csv')
edited = combined_feat.drop(columns=['ID','figo']) # 'shape_VoxelVolume', 'figo', 'ID'
cohort_norm=False
stage=False

# assign staging
values = {'1B':0, '1B1':0, '2A':0, '2B':0, '3A':1, '3B':1, '3C':1, '3C1':1, '3C2':1, '4A':1, '4B':1}
a = list(edited['figo2018'])
print(a, len(a))
b = []

for val in a:
    b.append(values[val])
edited['figo2018'] = b

# idify cohort to normalize by it...
vals = {'Neitherlands':0, 'TCGA':1}
a = list(edited ['cohort'])
b = []

for val in a:
    b.append(vals[val])

edited['cohort'] = b
edited = edited.astype('double')

# divide into targets/features...
# Do not normalize volume...

train_set, train_set_target = edited.drop(columns=['event','time'], axis=1), edited[['event','time','figo2018']]

a = list(train_set)
# run to drop shape features...
for lab in a:
    if 'shape' in str(lab):
        if 'VoxelVolume' not in str(lab):
            train_set=train_set.drop(columns=[str(lab)])
            print(f"{lab} 'Dropped'.")


a = list(train_set)

max_iter = 10
results = {}

for k in range(max_iter):

  try:
      train_set_, test_set_, train_set_target_,  test_set_target_  = train_test_split(train_set, train_set_target, test_size=0.2)
      train_set_, test_set_ = normalize(train_set_, test_set_)
      lab = get_signatures(train_set_, test_set_, train_set_target_,  test_set_target_)

  except Exception:
      train_set_, test_set_, train_set_target_,  test_set_target_  = train_test_split(train_set, train_set_target, test_size=0.2)
      train_set_, test_set_ = normalize(train_set_, test_set_)
      lab = get_signatures(train_set_, test_set_, train_set_target_,  test_set_target_)

  print(lab.keys())
  results[k] = lab
  print(f'Done {k}')

# save results...
ts = time.time()
date = datetime.datetime.fromtimestamp(ts).strftime("%Y_%m_%ds_%H%M%S")
save_obj(results,f'/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/-RADREGNORM_{date}_{max_iter}.pkl')
