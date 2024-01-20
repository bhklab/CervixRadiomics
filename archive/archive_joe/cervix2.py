import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pymrmre
import seaborn as sns
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
import pickle, json, time, datetime

def save_obj_json(obj, name):
    with open(name, 'w') as fp:
        json.dump(obj, fp)

def load_obj(name):
    with open(name, 'rb') as f:
        return pickle.load(f)

def save_obj(obj, name):
    with open(name, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

select = pd.read_csv('/Users/josephmarsilla/Cervix/Data/Features_with_all_ICC_above_0.5.csv')
combined_feat = pd.read_csv('/Users/josephmarsilla/Cervix/Combinied_data_for_Joe_2020_06_16.csv')

edited = combined_feat.drop(columns=['ID','figo']) # 'shape_VoxelVolume', 'figo', 'ID'

# Volume ONLY
# Volume + Radiomics
# Radiomics ONLY
# Stage ONLY
# Stage + Radiomics
# clinically they look at stage...
# take home message is we should conisder volume...
# radomics will add...

# assess for staging at the end...
# to test by stage...
values = {'1B':0, '1B1':0, '2A':0, '2B':0, '3A':1, '3B':1, '3C':1, '3C1':1, '3C2':1, '4A':1, '4B':1}
a = list(edited['figo2018'])
print(a, len(a))
b = []
for val in a:
    b.append(values[val])
edited['figo2018'] = b

# idify cohort to normalize by it...
vals = {'Neitherlands':1, 'TCGA':2}
a = list(edited['cohort'])
b = []
for val in a:
    b.append(vals[val])
edited['cohort'] = b

edited = edited.astype('double')

# list_ = []
# for i in range(len(edited.iloc[0])):
#     list_.append(0)

train_set, train_set_target = edited.drop(columns=['event','time'], axis=1), edited[['event','time','figo2018']]

results = {}

max_iter = 100

for itter_num in range(max_iter):

    train_set_, test_set_, train_set_target_,  test_set_target_  = train_test_split(train_set, train_set_target, test_size=0.2)

    counts = [1]
    result= { 1:[]}

    # stats = []
    # z-score normalization (by cohort)
    a = list(train_set_)
    ignore = ['cohort', 'figo2018']
    normalized_trainset = {}
    normalized_testset = {}

    for i, val in enumerate(a):
        b = []
        c = []
        if val not in ignore:

            for j in range(1,3):
                # normalize by cohort...
                me = train_set_[train_set_['cohort']==j][val].mean()
                stds = train_set_[train_set_['cohort']==j][val].std()
                b += list((train_set_[train_set_['cohort']==j][val] - me)/stds)
                c += list((test_set_[test_set_['cohort']==j][val] - me)/stds)

            normalized_trainset[val] = b
            normalized_testset[val] = c

    normalized_trainset['figo2018'] = train_set_['figo2018']
    normalized_testset['figo2018'] = test_set_['figo2018']
    train_set_ = pd.DataFrame.from_dict(normalized_trainset)
    test_set_ = pd.DataFrame.from_dict(normalized_testset)
    train_set_ = train_set_.astype('double')
    test_set_ = test_set_.astype('double')

    # stats.append((me,stds))
    # normalized trainset to dataframe...
    # include cohort type in the model...
    # Netherlands have the best survival rate...
    # TCGA had lowest survival rate because of stage...

    for j, c in enumerate(counts):
        cindices = []
        print("The current solution count is", c)
        for l in range(1, 41):

            # print("The current features selected length is", length)
            # print("Now we are testing the cox model with fixed")
            avg_c1, avg_c2 = 0.0, 0.0
            print('Hey')
            s1 = pymrmre.mrmr_ensemble(features=train_set_.drop(columns=['figo2018']), targets = train_set_target_.drop(columns=['time','figo2018']),
                                       solution_length=l, solution_count=c, fixed_features=['shape_VoxelVolume']) # 'shape_VoxelVolume', 'figo'


            for i in range(c):

                train_mrmr = train_set_[s1.iloc[0][i]]
                # what we'll be training on...
                train_mrmr['event'] = train_set_target_['event']
                train_mrmr['time'] = train_set_target_['time']
                # break up test set by stage...
                test_mrmr = test_set_[s1.iloc[0][i]]
                test_mrmr['event'] = test_set_target_['event']
                test_mrmr['time'] = test_set_target_['time']

                # can normalize it...
                cph1 = CoxPHFitter(penalizer=0.01) #penalizer=0.01
                cph1.fit(train_mrmr, duration_col='time', event_col='event', show_progress=False, step_size=0.25)
                # lower step size to a quarter ...
                # do this for different test cohorts...
                #try:
                c1 = concordance_index(test_mrmr['time'], -cph1.predict_partial_hazard(test_mrmr).values, test_mrmr['event'])

                avg_c1 += c1
                # avg_c2 += c2


            avg_c1 /= c
            # avg_c2 /= c


            cindices.append(avg_c1)

        result[c] = cindices
        print("The analysis of this solution count is all done!")
        print(cindices)
        # break

    # results.append(result)
    results[itter_num] = result
    print("All done!")
    # break

# avg_ci = 0
# for ins in range(len(results)):
#     avg_ci += np.array(results[ins][1])
#
# print(avg_ci/max_iter)

# save results...
ts = time.time()
date = datetime.datetime.fromtimestamp(ts).strftime("%Y_%m_%d_%H%M%S")
save_obj(results,f'/Users/josephmarsilla/Cervix/Data/VOLUMECOHORTNORM2_{date}_{max_iter}.pkl')
