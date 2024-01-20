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
cohort_norm=False
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
a = list(train_set)
# for lab in a:
#     if 'shape' in str(lab):
#         if 'VoxelVolume' not in str(lab):
#             train_set=train_set.drop(columns=[str(lab)])
#             print(f"{lab} 'Dropped'.")

results = {}
max_iter = 100
# label distribution ...
label = {}
countss = 0

for itter_num in range(max_iter):

    train_set_, test_set_, train_set_target_,  test_set_target_  = train_test_split(train_set, train_set_target, test_size=0.2)

    # try:
    #     assert len(test_set_[test_set_['figo2018'] == 0]) > 7
    #     assert len(test_set_[test_set_['figo2018'] == 1]) > 7
    #
    # except Exception:
    #
    #     print('Ensuring at least 5 of each stage in test set...')
    #     train_set_, test_set_, train_set_target_,  test_set_target_  = train_test_split(train_set, train_set_target, test_size=0.2)
    #
    #     a = len(test_set_[test_set_['figo2018'] == 0])
    #     b = len(test_set_[test_set_['figo2018'] == 1])
    #
    #     while a < 7 and b < 7:
    #         train_set_, test_set_, train_set_target_,  test_set_target_  = train_test_split(train_set, train_set_target, test_size=0.2)
    #         a = len(test_set_[test_set_['figo2018'] == 0])
    #         b = len(test_set_[test_set_['figo2018'] == 1])
    #
    #     print('Done search...')

    #     test_set, test_set_target = data_opc2.drop(['event','time'], axis=1), opc2_target
    #     fixed_features = ['Age','Sex','ECOG PS','Smoking Hx','Drinking hx','T','N','Stage']
    #     category_features = ['Sex','ECOG PS','Smoking Hx','Drinking hx','T','N','Stage']

    counts = [1]
    result= { 1:[]}

    # stats = []
    # z-score normalization (by cohort)
    a = list(train_set_)

    if cohort_norm:
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
        for l in range(2, 7):

            #print("The current features selected length is", length)
            #print("Now we are testing the cox model with fixed")

            avg_c1, avg_c2 = 0.0, 0.0
            print('Hey')

            s1 = pymrmre.mrmr_ensemble(features=train_set_.drop(columns=['figo2018']), targets = train_set_target_.drop(columns=['time','figo2018']),
                                       solution_length=l, solution_count=c, fixed_features=['cohort','shape_VoxelVolume'])

            # 'figo' # 'shape_VoxelVolume', keep cohort in there!! , 'shape_VoxelVolume' .drop(columns=[])

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

                # 2 year disease free survival, binarize
                # train_mrmr[train_mrmr['time'] < 2] = 0
                # train_mrmr[train_mrmr['time'] > 2] = 1
                # test_mrmr[test_mrmr['time'] < 2] = 0
                # test_mrmr[test_mrmr['time'] > 2] = 1

                # if stage:
                #     test_mrmr1 = test_mrmr[test_mrmr['figo2018']==0]
                #     test_mrmr2 = test_mrmr[test_mrmr['figo2018']==1]
                # else:
                #     # break up test set by stage...
                #     test_mrmr1 = test_set_[test_set_['figo2018']==0][s1.iloc[0][i]]
                #     test_mrmr2 = test_set_[test_set_['figo2018']==1][s1.iloc[0][i]]
                    # test_mrmr = test_set_[s1.iloc[0][i]]
    #                 for i, val in enumerate(list(s1.iloc[0][i])):
    #                     test_mrmr[val] = (test_mrmr[val] - vals[i][0])/vals[i][1]
                # test_mrmr2['event'] = test_set_target_[test_set_target_['figo2018']==1]['event']
                # test_mrmr2['time'] = test_set_target_[test_set_target_['figo2018']==1]['time']

                # print(test_mrmr.isnull().sum())

                # what if we run on pre-normalized data?
                # can normalize it...
                cph1 = CoxPHFitter() #penalizer=0.01
                cph1.fit(train_mrmr, duration_col='time', event_col='event', show_progress=False) # step_size=0.25
                print(cph1.concordance_index_)

                # if l == 2:
                #     # cphf.fit(train_mrmr[['event','time', 'shape_VoxelVolume']], duration_col='time', event_col='event', show_progress=False) # step_size=0.1
                #     print(cph1.concordance_index_)
                #     test_signature = -cph1.predict_partial_hazard(test_mrmr).values

                # export coeficcients of the fit ...
                coef = cph1.params_
                # print(coef['cohort'], coef)
                # ['coef']
                # get the median for each signature
                a = list(train_mrmr)
                exclude = ['cohort', 'time', 'event','shape_VoxelVolume']
                for c, val in enumerate(a):
                    # if val == 'shape_VoxelVolume':
                    #     train_mrmr[val] *= coef[val]
                    # if val == 'cohort':
                    #     train_mrmr[val] *= coef[val]
                    # leave cohort out
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

                # train_mrmr['Signature'] = (train_mrmr['Signature'] - train_mrmr['Signature'].mean()) /train_mrmr['Signature'].std()
                # test_mrmr['Signature'] = (test_mrmr['Signature']- train_mrmr['Signature'].mean()) /train_mrmr['Signature'].std()
                print(test_mrmr['Signature'])
                # fit a cox model with that signature
                # run it on test...
                cphf =  CoxPHFitter() #penalizer=0.01
                if l == 2:
                    cphf.fit(train_mrmr[['event','time', 'cohort', 'shape_VoxelVolume']], duration_col='time', event_col='event', show_progress=False) # step_size=0.1
                    print(cphf.concordance_index_)
                    test_signature = -cphf.predict_partial_hazard(test_mrmr[['event','time', 'cohort', 'shape_VoxelVolume']]).values
                if l>2:
                    cphf.fit(train_mrmr[['event','time', 'shape_VoxelVolume', 'Signature']], duration_col='time', event_col='event', show_progress=False) # step_size=0.1
                    print(cphf.concordance_index_)
                    test_signature = -cphf.predict_partial_hazard(test_mrmr[['event','time', 'shape_VoxelVolume', 'Signature']]).values


                # test_signature2 = -cphf.predict_survival_function(test_mrmr[['event','time','Signature']]).values
                print(test_signature)
                # break
                # calculate signature for test set...
                # then run with the above mode

                # take the median after we add together...
                # training to get signature ...
                # use that signature and cut point in the test (dicotomized = divide by median training signature)

                # train fit test coefficient
                # test fit to test signatures
                # run the whole thing with one iteration ...

                # train_median = train_mrmr.median()
                # print(train_median)

                # divide test_mrmr by median?? to get 'dichotomized signature'
                # get signature for all patients
                # divide by median ...
                # for each patient we get a signature score ...
                # in training, per iteration...
                # different signature with different median ...
                # we're seeing if the method works, what is the optimal number of features?
                # one signature then to apply on pmh ...
                # dichotomized_mrmr = test_mrmr.copy()
                # for c, val in enumerate(a):
                #     if val not in exclude:
                #         dichotomized_mrmr[val] /= train_median[val]
                # dico_signature = -cph1.predict_partial_hazard(dichotomized_mrmr).values

                # lower step size to a quarter ...
                # do this for different test cohorts...
                # try:

                c1 = concordance_index(test_mrmr['time'], test_signature, test_mrmr['event'])
                # c2 = concordance_index(test_mrmr['time'], dico_signature, test_mrmr['event'])

                #
                # except Exception:
                #     print('Bad run c1')
                #     c1 = 0
                # try:
                #     c2 = concordance_index(test_mrmr2['time'], -cph1.predict_partial_hazard(test_mrmr2).values, test_mrmr2['event'])
                # except Exception:
                #     print('Bad run c2')
                #     c2 = 0

                avg_c1 += c1
                # avg_c2 += c2

                #print("The current c-index for model-1 is", c1)

            avg_c1 /= c
            # avg_c2 /= c


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
save_obj(results,f'/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/VOLRADREGNORM_{date}_{max_iter}.pkl')
save_obj(label,f'/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/VOLRADLABELS_{date}_{max_iter}.pkl')
