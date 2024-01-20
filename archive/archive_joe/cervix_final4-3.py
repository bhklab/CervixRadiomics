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
import pickle, json, time, datetime, argparse, warnings
# from google.colab import drive
# drive.mount('/content/drive')
# %matplotlib inline
# set working directory
# %cd drive/My\ Drive/notebooks/Cervix/
# %cd /Users/josephmarsilla/Google\ Drive/notebooks/Cervix/
def get_args():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--max-iter', default='10', type=int,
                        help='Max Iterations')
    parser.add_argument('--n-solutions', default='100', type=int,
                        help='Number of solutions for ensemble.')
    parser.add_argument('--choose-data', default='both', type=str, help='both is default. (TCGA == 1) (Norway == 0)')
    # parser.add_argument('--sum', dest='accumulate', action='store_const',
    #                     const=sum, default=max,
    #                     help='sum the integers (default: find the max)')

    args = parser.parse_args()
    return args

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
def get_signatures(train_set_, test_set_, train_set_target_,  test_set_target_, solutions=100, data='both'):
      lab = {}
      if data == 'both':
        start = 2
        end = 7
      else:
        start = 1
        end = 6

      for l in range(start,end):
        doc = {}
        avg_c1, avg_c2 = 0.0, 0.0

        print(f'Starting {l}') # .drop(columns=['shape_VoxelVolume'])

        if data == 'both':
            drop_cols = ['shape_VoxelVolume', 'figo2018']
            fixed = ['cohort']
        else:
            drop_cols = ['shape_VoxelVolume', 'figo2018', 'cohort']
            fixed = []

        doc = {}
        avg_c1, avg_c2 = 0.0, 0.0
        s = pymrmre.mrmr_ensemble(features=train_set_.drop(columns=drop_cols), targets = train_set_target_.drop(columns=['time','figo2018']),
                                    solution_length=l, solution_count=1, fixed_features=fixed)

        s1 = pymrmre.mrmr_ensemble(features=train_set_.drop(columns=drop_cols), targets = train_set_target_.drop(columns=['time','figo2018']),
                                    solution_length=l, solution_count=solutions, fixed_features=fixed)

        s2 = pymrmre.mrmr_ensemble(features=train_set_.drop(columns=drop_cols), targets = train_set_target_.drop(columns=['time','figo2018']),
                                    solution_length=l, solution_count=solutions*5, fixed_features=fixed)

        feature_list0 = list(s.iloc[0])
        feature_list = list(s1.iloc[0])
        feature_list2 = list(s2.iloc[0])
        doc['feature_list_1'] = feature_list
        doc['feature_list_100'] = feature_list
        doc['feature_list_500'] = feature_list2
        doc['num_features'] = l - 2 # in addition to volume...

        # for every solution in the ensemble calculate a train/test ci...
        # should be length solution_count
        test_ci = []
        train_ci = []
        # coeficients used in enseble...
        coefs = []

        train_mrmr = train_set_.copy()
        train_mrmr['event'] = train_set_target_['event']
        train_mrmr['time'] = train_set_target_['time']
        test_mrmr = test_set_.copy()
        test_mrmr['event'] = test_set_target_['event']
        test_mrmr['time'] = test_set_target_['time']

        # set_features = ['cohort']
        if data == 'both':
            penalizer = None
        else:
            penalizer = 0.01

        for z, set_ in enumerate(feature_list0):

            exclude = ['cohort', 'time', 'event','shape_VoxelVolume', 'figo2018']

            # fit cox model with chosen features
            cph1 = CoxPHFitter(penalizer=penalizer) #penalizer=0.01
            cph1.fit(train_mrmr[set_ + ['time', 'event']], duration_col='time', event_col='event', show_progress=False) # step_size=0.25
            # doc['train_ci'] = cph1.concordance_index_
            coef = cph1.params_
            # coefs.append(coef)
            # doc['cox_coef'] = coef

    #             a = list(b)
    #             choice = ['shape_VoxelVolume']

            for c, val in enumerate(set_):

              # if val in choice:
              #     b[val] *= coef[val]
              #     btest[val] *= coef[val]

              if val not in exclude:
                  v =  train_mrmr[val]*coef[val]
                  v_test =  test_mrmr[val]*coef[val]
                  # calculates signature from all the features ...
                  try:
                      train_mrmr['Signature0'] += v
                      test_mrmr['Signature0'] += v_test
                  except Exception:
                      train_mrmr['Signature0'] = v
                      test_mrmr['Signature0'] = v_test

            if z%10 == 0:
                print(f'Done solution {z} for {l} features.')

        print('Completed 1 solutions. Starting 100.')

        for z, set_ in enumerate(feature_list):

            exclude = ['cohort', 'time', 'event','shape_VoxelVolume', 'figo2018']

            # fit cox model with chosen features
            cph1 = CoxPHFitter(penalizer=penalizer) #penalizer=0.01
            cph1.fit(train_mrmr[set_ + ['time', 'event']], duration_col='time', event_col='event', show_progress=False) # step_size=0.25
            # doc['train_ci'] = cph1.concordance_index_
            coef = cph1.params_

            for c, val in enumerate(set_):

              # if val in choice:
              #     b[val] *= coef[val]
              #     btest[val] *= coef[val]

              if val not in exclude:
                  v =  train_mrmr[val]*coef[val]
                  v_test =  test_mrmr[val]*coef[val]
                  # calculates signature from all the features ...
                  try:
                      train_mrmr['Signature'] += v
                      test_mrmr['Signature'] += v_test
                  except Exception:
                      train_mrmr['Signature'] = v
                      test_mrmr['Signature'] = v_test

            if z%10 == 0:
                print(f'Done solution {z} for {l} features.')

        print('Completed 100 solutions. Starting 500.')

        for z, set_ in enumerate(feature_list2):

            exclude = ['cohort', 'time', 'event','shape_VoxelVolume', 'figo2018']

            # fit cox model with chosen features
            cph1 = CoxPHFitter(penalizer=penalizer) #penalizer=0.01
            cph1.fit(train_mrmr[set_ + ['time', 'event']], duration_col='time', event_col='event', show_progress=False) # step_size=0.25
            # doc['train_ci'] = cph1.concordance_index_
            coef = cph1.params_
            # coefs.append(coef)
            # doc['cox_coef'] = coef

    #             a = list(b)
    #             choice = ['shape_VoxelVolume']

            for c, val in enumerate(set_):

              # if val in choice:
              #     b[val] *= coef[val]
              #     btest[val] *= coef[val]

              if val not in exclude:
                  v =  train_mrmr[val]*coef[val]
                  v_test =  test_mrmr[val]*coef[val]
                  # calculates signature from all the features ...
                  try:
                      train_mrmr['Signature2'] += v
                      test_mrmr['Signature2'] += v_test
                  except Exception:
                      train_mrmr['Signature2'] = v
                      test_mrmr['Signature2'] = v_test

            if z%10 == 0:
                print(f'Done solution {z} for {l} features.')

        if data=='both':
            chosen0 = ['event','time', 'cohort','Signature0']
            chosen = ['event','time', 'cohort', 'Signature']
            chosen2 = ['event','time', 'cohort','Signature2']
        else:
            chosen0 = ['event','time', 'Signature0']
            chosen = ['event','time', 'Signature']
            chosen2 = ['event','time', 'Signature2']


        ##################################
        b = train_mrmr[chosen0]
        btest = test_mrmr[chosen0]

        # divide by total ensemble number
        b['Signature0'] /= len(feature_list0)
        btest['Signature0'] /= len(feature_list0)
        ##################################
        # fit final model for ensemble...
        cphf =  CoxPHFitter() # penalizer=0.01
        cphf.fit(b, duration_col='time', event_col='event', show_progress=False) # step_size=0.1
        # calculate train concordance
        c1 = cphf.concordance_index_
        train_ci.append(c1)

        # use signature to calculate test concordance...
        test_signature = -cphf.predict_partial_hazard(btest).values
        c2 = concordance_index(btest['time'], test_signature, btest['event'])
        test_ci.append(c2)

        doc['train_ci_1'] = c1
        doc['test_ci_1'] = c2

        print(f'Done 1 run {l}. Average Train CI == {c1}. Average Test CI == {c2}')

        ##################################
        b = train_mrmr[chosen]
        btest = test_mrmr[chosen]

        # divide by total ensemble number
        b['Signature'] /= len(feature_list)
        btest['Signature'] /= len(feature_list)
        ##################################
        # fit final model for ensemble...
        cphf =  CoxPHFitter() # penalizer=0.01
        cphf.fit(b, duration_col='time', event_col='event', show_progress=False) # step_size=0.1
        # calculate train concordance
        c1 = cphf.concordance_index_
        train_ci.append(c1)

        # use signature to calculate test concordance...
        test_signature = -cphf.predict_partial_hazard(btest).values
        c2 = concordance_index(btest['time'], test_signature, btest['event'])
        test_ci.append(c2)

        doc['train_ci_100'] = c1
        doc['test_ci_100'] = c2

        print(f'Done 100 run {l}. Average Train CI == {c1}. Average Test CI == {c2}')

        ########################
        b = train_mrmr[chosen2]
        btest = test_mrmr[chosen2]

        # divide by total ensemble number
        b['Signature2'] /= len(feature_list2)
        btest['Signature2'] /= len(feature_list2)
        ########################

        # fit final model for ensemble...
        cphf =  CoxPHFitter() # penalizer=0.01
        cphf.fit(b, duration_col='time', event_col='event', show_progress=False) # step_size=0.1
        # calculate train concordance
        c1 = cphf.concordance_index_
        train_ci.append(c1)

        # use signature to calculate test concordance...
        test_signature = -cphf.predict_partial_hazard(btest).values
        c2 = concordance_index(btest['time'], test_signature, btest['event'])
        test_ci.append(c2)

        doc['train_ci_500'] = c1
        doc['test_ci_500'] = c2
        # doc['cox_coefs'] = coefs

        lab[l-1] = doc
        print(f'Done 500 run {l}. Average Train CI == {c1}. Average Test CI == {c2}')

      return lab

args = get_args()

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
d = args.choose_data
if d != 'both':
    # only choose Netherlands
    warnings.warn('Using Norway Data ONLY')
    edited = edited[edited['cohort'] == 0]

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

max_iter = args.max_iter
results = {}

for k in range(max_iter):

  try:
      train_set_, test_set_, train_set_target_,  test_set_target_  = train_test_split(train_set, train_set_target, test_size=0.2)
      train_set_, test_set_ = normalize(train_set_, test_set_)
      lab = get_signatures(train_set_, test_set_, train_set_target_,  test_set_target_, max_iter,d)

  except Exception:
      train_set_, test_set_, train_set_target_,  test_set_target_  = train_test_split(train_set, train_set_target, test_size=0.2)
      train_set_, test_set_ = normalize(train_set_, test_set_)
      lab = get_signatures(train_set_, test_set_, train_set_target_,  test_set_target_, max_iter,d)

  print(lab.keys())
  results[k] = lab
  print(f'Done {k}')

# save results...
ts = time.time()
date = datetime.datetime.fromtimestamp(ts).strftime("%Y_%m_%ds_%H%M%S")
save_obj(results,f'/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/-RADREGNORM_{date}_{max_iter}_{d}.pkl')
