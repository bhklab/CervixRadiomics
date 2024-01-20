import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pymrmre
import seaborn as sns
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from sklearn.model_selection import train_test_split

combined_feat = pd.read_csv('/Users/josephmarsilla/Cervix/Data/final_train_edited.csv')
feats = list(combined_feat)[3:]
edited = combined_feat.drop(columns=['PatientID'])
edited = edited.astype('double')
list_ = []
for i in range(len(edited.iloc[0])):
    list_.append(0)
train_set, train_set_target = edited.drop(['event','time'], axis=1), edited[['event','time']]

results = []

for itter_num in range(5):
    train_set, test_set, train_set_target,  test_set_target  = train_test_split(train_set, train_set_target, test_size=0.2)
    #test_set, test_set_target = data_opc2.drop(['event','time'], axis=1), opc2_target
    #     fixed_features = ['Age','Sex','ECOG PS','Smoking Hx','Drinking hx','T','N','Stage']
    #     category_features = ['Sex','ECOG PS','Smoking Hx','Drinking hx','T','N','Stage']
    counts = [1]
    result= { 1:[[], [],[],[]]}

    for j, c in enumerate(counts):
        cindices = [[], [],[],[]]
        print("The current solution count is", c)
        for l in range(10, 41):

            #print("The current features selected length is", length)
            #print("Now we are testing the cox model with fixed")

            avg_c1, avg_c2 = 0.0, 0.0
            print('Hey')

            s1 = pymrmre.mrmr_ensemble(features=train_set, targets = train_set_target,
                                       category_features = feats,
                                       solution_length=l, solution_count=c)

#             s1 = mrmr.mrmr_ensemble_survival(features = train_set,
#                                              targets = train_set_target,
#                                              category_features = category_features,
#                                              fixed_features = fixed_features,
#                                              solution_count = c,
#                                              solution_length = l)


            for i in range(c):

                train_mrmr = train_set[s1.iloc[0][i]]
                train_mrmr['event'] = train_set_target['event']
                train_mrmr['time'] = train_set_target['time']
                test_mrmr = test_set[s1.iloc[0][i]]
                test_mrmr['event'] = test_set_target['event']
                test_mrmr['time'] = test_set_target['time']

                cph1 = CoxPHFitter()
                cph1.fit(train_mrmr, duration_col='time', event_col='event', show_progress=False)
                c1 = concordance_index(test_mrmr['time'], -cph1.predict_partial_hazard(test_mrmr).values, test_mrmr['event'])
                avg_c1 += c1
                #print("The current c-index for model-1 is", c1)

            avg_c1 /= c


            cindices[0].append(avg_c1)

        result[c] = cindices
        print("The analysis of this solution count is all done!")
        print(cindices)

    results.append(result)
    print("All done!")
