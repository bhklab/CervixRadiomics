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

select = pd.read_csv('/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/Features_with_all_ICC_above_0.5.csv')
combined_feat = pd.read_csv('/Users/josephmarsilla/Google Drive/notebooks/Cervix/Combinied_data_for_Joe_2020_06_16.csv')
edited = combined_feat.drop(columns=['ID','figo']) # 'shape_VoxelVolume', 'figo', 'ID'
cohort_norm=True
stage=False

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
# they would like to include cervices in the model

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

# list_ = []
# for i in range(len(edited.iloc[0])):
#     list_.append(0)

train_set, train_set_target = edited.drop(columns=['event','time'], axis=1), edited[['event','time','figo2018']]
# drop everything correlated to shape...
# need to find BEST radiomic features
# include and disclude ...
a = list(train_set)
for lab in a:
    if 'shape' in str(lab):
        train_set=train_set.drop(columns=[str(lab)])
        print(f"{lab} 'Dropped'.")

results = {}
max_iter = 1000
# label distribution ...
label = {}
countss = 0

for itter_num in range(max_iter):

    train_set_, test_set_, train_set_target_,  test_set_target_  = train_test_split(train_set, train_set_target, test_size=0.2)

    try:
        assert len(test_set_[test_set_['figo2018'] == 0]) > 7
        assert len(test_set_[test_set_['figo2018'] == 1]) > 7

    except Exception:

        print('Ensuring at least 5 of each stage in test set...')
        train_set_, test_set_, train_set_target_,  test_set_target_  = train_test_split(train_set, train_set_target, test_size=0.2)

        a = len(test_set_[test_set_['figo2018'] == 0])
        b = len(test_set_[test_set_['figo2018'] == 1])

        while a < 7 and b < 7:
            train_set_, test_set_, train_set_target_,  test_set_target_  = train_test_split(train_set, train_set_target, test_size=0.2)
            a = len(test_set_[test_set_['figo2018'] == 0])
            b = len(test_set_[test_set_['figo2018'] == 1])

        print('Done search...')

    #     test_set, test_set_target = data_opc2.drop(['event','time'], axis=1), opc2_target
    #     fixed_features = ['Age','Sex','ECOG PS','Smoking Hx','Drinking hx','T','N','Stage']
    #     category_features = ['Sex','ECOG PS','Smoking Hx','Drinking hx','T','N','Stage']

    counts = [1]
    result= { 1:[]}

    # stats = []
    # z-score normalization (by cohort)
    a = list(train_set_)

    if cohort_norm:
        a = list(train_set_)
        ignore = ['cohort', 'figo2018']
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
        ignore = [] # 'cohort', 'figo2018'
        for i, val in enumerate(a):
            if val not in ignore:
                me = train_set_[val].mean()
                stds = train_set_[val].std()
                train_set_[val] = (train_set_[val] - me)/stds
                test_set_[val]  = (test_set_[val] - me)/stds

    # stats.append((me,stds))
    # normalized trainset to dataframe...
    # include cohort type in the model...
    # Netherlands have the best survival rate...
    # TCGA had lowest survival rate because of stage...
    lab = {}
    for j, c in enumerate(counts):
        cindices = []
        print("The current solution count is", c)
        for l in range(2, 10):

            #print("The current features selected length is", length)
            #print("Now we are testing the cox model with fixed")

            avg_c1, avg_c2 = 0.0, 0.0
            print('Hey')

            s1 = pymrmre.mrmr_ensemble(features=train_set_.drop(columns=['figo2018']), targets = train_set_target_.drop(columns=['time','figo2018']),
                                       solution_length=l, solution_count=c, fixed_features=['cohort'])

            # 'figo' # 'shape_VoxelVolume', keep cohort in there!!

            lab[l-1] = list(s1.iloc[0])[0][1:] # assuming first column is

#             s1 = mrmr.mrmr_ensemble_survival(features = train_set,
#                                              targets = train_set_target,
#                                              category_features = category_features,
#                                              fixed_features = fixed_features,
#                                              solution_count = c,
#                                              solution_length = l)

#             x = train_set_.values #returns a numpy array
#             min_max_scaler = preprocessing.MinMaxScaler()
#             x_scaled = min_max_scaler.fit_transform(x)
#             train_set_ = pd.DataFrame(x_scaled, columns=train_set_.columns)
#             test_scaled = min_max_scaler.transform(test_set_.values)
#             test_set_ = pd.DataFrame(test_scaled, columns=train_set_.columns)


            for i in range(c):

                train_mrmr = train_set_[s1.iloc[0][i]]
                # what we'll be training on...
                train_mrmr['event'] = train_set_target_['event']
                train_mrmr['time'] = train_set_target_['time']

                test_mrmr = test_set_[s1.iloc[0][i]]
                test_mrmr['event'] = test_set_target_['event']
                test_mrmr['time'] = test_set_target_['time']
                # what if we run on pre-normalized data?
                # can normalize it...
                cph1 = CoxPHFitter(penalizer=0.01) #penalizer=0.01
                cph1.fit(train_mrmr, duration_col='time', event_col='event', show_progress=False, step_size=0.25)
                # lower step size to a quarter ...
                # do this for different test cohorts...
                coef = cph1.params_
                # get the median for each signature
                a = list(train_mrmr)
                exclude = ['cohort', 'time', 'event']
                for c, val in enumerate(a):
                    # leave cohort out
                    if val not in exclude:
                        train_mrmr[val] *=  coef[val]

                train_median = train_mrmr.median()
                print(train_median)

                # divide test_mrmr by median?? to get 'dichotomized signature'
                dichotomized_mrmr = test_mrmr.copy()
                for c, val in enumerate(a):
                    if val not in exclude:
                        dichotomized_mrmr[val] /= train_median[val]
                # try:
                c1 = concordance_index(dichotomized_mrmr['time'], -cph1.predict_partial_hazard(dichotomized_mrmr).values, dichotomized_mrmr['event'])
                c2 = 0

                avg_c1 += c1
                avg_c2 += c2

                #print("The current c-index for model-1 is", c1)

            avg_c1 /= c
            avg_c2 /= c


            # cindices.append([avg_c1, avg_c2])
            cindices.append(avg_c1)

        result[c] = cindices
        print("The analysis of this solution count is all done!")
        print(cindices)
        # break

    # results.append(result)
    results[itter_num] = result
    label[itter_num] = lab
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
save_obj(results,f'/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/OLD_DICRADREGNORM_{date}_{max_iter}.pkl')
save_obj(label,f'/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/OLD_DICLABELS_{date}_{max_iter}.pkl')
