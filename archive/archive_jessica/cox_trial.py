import numpy as np
import pickle, os, datetime, time, pathlib
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines import CoxPHFitter
import pymrmre
%matplotlib inline
%config InlineBackend.figure_format = 'retina'
from matplotlib import pyplot as plt
# use combined features to choose data
help(pymrmre)# load from a pandas dataframe
combined_features = pd.read_csv('combined_features.csv')
selected = list(combined_features.head())[1:]
list_ = []

edited = combined_features.drop(columns=['PatientID'])
edited = edited.astype('double')
for i in range(len(edited.iloc[0])):
    list_.append(0)

# how do we set fixed feature count?
# solution_count=10
# solution_length=10
# fixed_features=3

solns = {}
# to keep power...
for i in range(6,16):

    for j in range(0,6):
        solutions = pymrmre.mrmr_ensemble(features=edited, feature_types=list_,
                                          target_features=selected, solution_length=i,
                                          solution_count=10, fixed_feature_count=j)
        solns[f'{i},{j}'] = solutions

    print(i)

output = open('solutions_cervix_march_32020.pkl', 'wb')
pickle.dump(solns, output)

# complete the same analysis with normalized version of features...
solns['15,5'][0][0]

# this normalizes columns by subtracting std and dividing mean...
import sklearn
from sklearn import preprocessing as pre
x = edited.values #returns a numpy array
scaler = pre.StandardScaler()
x_scaled = scaler.fit_transform(x)
edited_norm = pd.DataFrame(x_scaled, columns=edited.columns)
# check head ...
edited_norm.head()

# import actual training data, normalize
final = pd.read_csv('/Users/jmarsill/Desktop/AI/Cervix/combined_features_2.csv')
final = final.drop('PatientID', axis=1)
final_ = final.copy()
y = final_.values
y_scaled = scaler.fit_transform(y)
final_norm = pd.DataFrame(y_scaled, columns=final_.columns)
# check head
final_norm.head()
final_norm['dftime_rad'] = final['dftime_rad']
final_norm['dfcens'] = final['dfcens']
final_norm.head()

# redue for normalized features
aa = final_norm.copy()
aa = aa.drop('dftime_rad', axis=1)
aa = aa.drop('dfcens', axis=1)
aa = aa.astype('double')
list_ = []
selected = list(aa.head())
for i in range(len(aa.iloc[0])):
    list_.append(0)
len(list_)
# check aa head...
aa.head()

norm_solns = {}
for i in range(6,16):
    for j in range(0,6):
        solutions = pymrmre.mrmr_ensemble(features=aa, feature_types=list_,
                                          target_features=selected, solution_length=i,
                                          solution_count=10, fixed_feature_count=j)
        norm_solns[f'{i},{j}'] = solutions
    print(i)
output = open('solutions_cervix_march_32020.pkl', 'wb')
pickle.dump(solns, output)

final_norm.head()
# run the exact

# redo model fam...
# select using MRMRE
# dont have to fix features
# all the Features
# 10, 20, 50, 100, 200
# apply cox models
# after training cox model,
# use cross validation...
# 80, 20 and 5 fold cross validation
# need to select the features from the training set
# select all the features from the training

# dictionary to save trained models to ...
dic_ = {}

#####################################
# run this on 10 models...
a = np.arange(0,758,1)
b = np.arange(0,10,1)

for i in range(100):
    a_c = np.random.choice(a)
    b_c = np.random.choice(b)
    choice = solutions[a_c][b_c]
    final_ = ['dftime_rad', 'dfcens']
    final_ += choice.copy()
    to_pick = final[final_]
    cph = CoxPHFitter(penalizer=0)
    cph.fit(to_pick, 'dftime_rad', 'dfcens', step_size=0.25)
    cph.print_summary()
    dic_[f'{a_c},{b_c}'] = cph

# '375,7' concordance of 0.78...
# hparams: 10,10,3
# '512,2'
# hparams: 15,10,5

len(list(dic_.keys()))
# save as pkl
import pickle
output = open('cph_15_5_115.pkl', 'wb')
pickle.dump(dic_, output)
dic

####################################
# running this on all models...
for i, soln in enumerate(solutions):
    for j, sol in enumerate(soln):
        choice = sol
        final_ = ['dftime_rad', 'dfcens']
        final_ += choice.copy()
        to_pick = final[final_]
        cph = CoxPHFitter(penalizer=0)
        cph.fit(to_pick, 'dftime_rad', 'dfcens', step_size=0.25)
        dic[f'{i},{j}'] = cph

    print(i)

    if i%10==0:
        cph.print_summary()
        break


# have to
# save ensemble...
# cph.plot_covariate_groups('TotalCharges', groups=[0,4000])


# final_.append('dftime_rad')
# final_.append('dfcens')
# final_
# solutions.head()
# features in the MRMR ensemble...
# solutions[0][5]
# %matplotlib inline
# %config InlineBackend.figure_format = 'retina'
# from matplotlib import pyplot as plt
# plt.style.use('bmh')


# to_pick
# print summaries of cox model
# dir(cph)
# cph.hazard_ratios_
# cph.baseline_survival_
# cph.plot()

# cph.predict_survival_function(to_pick).plot()

# cph.plot_covariate_groups('prio', [0, 2, 4, 6, 8, 10], cmap='coolwarm')


# # can run lifelines model with these features...
# # have to correlate feature data to event data...
#
# list(combined_features['PatientID'])
#
# # T = data["duration"]
# # E = data["observed"]
# #
# # kmf.fit(T, event_observed=E)
#
# netherlands_events = pd.read_csv('netherlands_events.csv')
# matching_neth  = pd.read_csv('/Users/jmarsill/Desktop/AI/Cervix/matching netherlands.csv')
# a = list(netherlands_events['PatientID'])
# a = [f'Patient{id}' for i, id in enumerate(a)]
# len(a)
# len(netherlands_events)
# neth_filt = matching_neth[matching_neth['Pt ID'].isin(a)].reset_index()['ID']
# len(neth_filt)
# neth_filt
#
#
# len(neth_filt)
# # netherlands_events
#
# tcga_events = pd.read_csv('tcga_events.csv')
# list(tcga_events[0])
# tcga_events['PatientID']
# # featues with events...
# len(list(tcga_events['PatientID']))
# tcga_events
# b = list(tcga_events['PatientID'])
# b
# len(b)
# TCGA = combined_features[combined_features['PatientID'].isin(b)]
# len(list(TCGA['PatientID']))
# filtered = [bs for bs in list(TCGA['PatientID']) if bs not in b]
# filtered
#
# # this puts everything out of order...
# # then we can run the cox model
# netherlands_events
# NETH = combined_features[combined_features['PatientID'].isin(neth_filt)].reset_index()
# # sort both events and extracted features by patientID order
# cf = NETH.set_index('PatientID')
# cf = cf.reindex(index=neth_filt)
# cf = cf.reset_index()
#
# ne = netherlands_events.set_index('PatientID')
# ne = ne.reindex(index=neth_filt)
# ne = ne.reset_index()
# netherlands_events
# neth_filt
#
# # now inserting cols should work...
# NETH.insert(2, 'dftime_rad', ne['dftime_rad'])
# NETH.insert(3, 'dfcens', ne['dfcens'])
# NETH
#
# TCGA = TCGA.reset_index()
# TCGA.drop(columns=['index'])
# TCGA.insert(2, 'dftime_rad', tcga_events['dftime_rad'])
# TCGA
# # tcga_events
#
# combined_events = pd.read_csv('combined_events.csv')
# # combined_events
#
#
# combined_features
# # combined_features
#
#
# # code allows us to make combined_features
# netherlands_features = pd.read_csv('neitherlands_radiomics_2019_09_10.csv')
# # len(netherlands_features)
#
# tcga_features = pd.read_csv('radiomics_features_GTV_MR_Diagnostic_TCGA.csv')
# tcga_features.head()
#
# # features with ICC abouve 0.5
# select_features = pd.read_csv('Features_with_all_ICC_above_0.5.csv')
# select = ['PatientID'] + select_features['Feature'].to_list()
# select
#
# # select reduced feature set
# neth_reduced = netherlands_features[select]
# select_ = [sel.replace('.', '-') for i, sel in enumerate(select) if 'log' or 'wavelet' in sel]
# select_
# tcga_reduced = tcga_features[select_]
# tcga_reduced.columns = select
#
# len(tcga_reduced)
# len(neth_reduced)
# len(list(tcga_reduced))
# len(list(neth_reduced))
# len(select)

# export data to csv, use mrmre package in R for feature selection
# then run a cox model baby...
# combined = pd.concat([tcga_reduced, neth_reduced])
# combined_events = pd.concat([tcga_events, netherlands_events])
# combined.to_csv('combined_features.csv', encoding='utf-8', index=False)
# combined_events.to_csv('combined_events.csv', encoding='utf-8', index=False)
