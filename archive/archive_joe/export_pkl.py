#######################
# summary of results...
#######################
import pandas as pd
import os, pickle
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

def load_obj(name):
    with open(name, 'rb') as f:
        return pickle.load(f)

# STAGE + VOL
# a = load_obj('/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/-VOLSTAGEREGNORM_2020_10_22s_124035_10.pkl')
# a = load_obj('/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/-VOLSTAGEREGNORM_2020_10_22s_092412_1000.pkl')
# a = load_obj('/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/-VOLSTAGEREGNORM_2020_09_01s_114231_100.pkl')
# a = load_obj('/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/NORONLY-VOLSTAGEREGNORM_2020_10_22s_130120_1000.pkl')
# a = load_obj('/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/NORONLY-VOLSTAGEREGNORM_2020_10_22s_131330_100.pkl')
a = load_obj('/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/NORONLY-VOLSTAGEREGNORM_2020_10_22s_131553_10.pkl')

# STAGE
# a = load_obj('/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/-STAGEREGNORM_2020_10_22s_122624_1000.pkl')
# a = load_obj('/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/-STAGEREGNORM_2020_09_01s_110614_100.pkl')
# a = load_obj('/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/-STAGEREGNORM_2020_10_22s_124030_10.pkl')
# a = load_obj('/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/NORONLY-STAGEREGNORM_2020_10_22s_131332_1000.pkl')
# a = load_obj('/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/NORONLY-STAGEREGNORM_2020_10_22s_131735_100.pkl')
# a = load_obj('/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/NORONLY-STAGEREGNORM_2020_10_22s_132406_10.pkl')

# VOL
# a = load_obj('/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/-VOLREGNORM_2020_09_01s_104935_100.pkl')
# a = load_obj('/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/-VOLREGNORM_2020_10_22s_122738_1000.pkl')
# a = load_obj('/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/-VOLREGNORM_2020_10_22s_124028_10.pkl')
# a = load_obj('/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/NORONLY-VOLREGNORM_2020_10_22s_124631_1000.pkl')
# a = load_obj('/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/NORONLY-VOLREGNORM_2020_10_22s_131320_100.pkl')
# a = load_obj('/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/NORONLY-VOLREGNORM_2020_10_22s_131537_10.pkl')

# RAD
# a = load_obj('/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/-RADREGNORM_2020_09_01s_111535_100.pkl')
# a = load_obj('/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/-RADREGNORM_2020_10_22s_124032_10.pkl')
# a = load_obj('/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/-RADREGNORM_2020_10_22s_122518_1000.pkl')
# a = load_obj('/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/NORONLY-RADREGNORM_2020_10_22s_132826_1000.pkl')
# a = load_obj('/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/NORONLY-RADREGNORM_2020_10_22s_133140_100.pkl')
a = load_obj('/Users/josephmarsilla/Google Drive/notebooks/Cervix/Data/NORONLY-RADREGNORM_2020_10_22s_135343_10.pkl')

iter = 10
file_name = 'norway_rad' + '_' + str(iter) + '_'
train_plots = [[], [], [], [], []]
test_plots =  [[], [], [], [], []]
dico_train_plots = [[], [], [], [], []]
dico_test_plots = [[], [], [], [], []]

# features = []

for i in range(len(a)):
    for j in range(len(a[0])):
        # add 1 to j for TCGA+NORWAY added cohort to model...
        test_plots[j].append(a[i][j]['test_ci'])
        train_plots[j].append(a[i][j]['train_ci'])
        dico_train_plots[j].append(a[i][j]['dicho_train_ci'])
        dico_test_plots[j].append(a[i][j]['dicho_test_ci'])

        # if j == 3:
        #     features += a[i][j+1]['features']

    print(f'Done {i}')

# names = ['', 'One Feature', 'Two features', 'Three Features', 'Four Features']
# name = 'Stage'

names = [ 'One Feature', 'Two features', 'Three Features', 'Four Features', 'Five Features']
name = 'Radiomics'

for i in range(len(train_plots)):

    plot = []
    chosen = ['Train', 'Test', 'DTrain', 'DTest']
    plot.append(train_plots[i])

    chosen[0] += f" ({str(np.round(np.mean(train_plots[i]),3))})"
    plot.append(test_plots[i])
    chosen[1] += f" ({str(np.round(np.mean(test_plots[i]),3))})"
    plot.append(dico_train_plots[i])
    chosen[2] += f" ({str(np.round(np.mean(dico_train_plots[i]),3))})"
    plot.append(dico_test_plots[i])
    chosen[3] += f" ({str(np.round(np.mean(dico_test_plots[i]),3))})"

    fig = plt.figure()
    plt.title(f"CI for Disease Free Survival ({name} + {names[i]})")
    plt.boxplot(plot)
    plt.xticks(range(1, len(chosen) + 1), chosen, rotation="horizontal")
    plt.ylabel("CI")
    plt.xlabel("Run Name")
    fig.savefig(f"{file_name}{names[i].replace(' ', '')}.png" )
    print(f'Done {names[i]}')

# feature plot, showing features selected by MRMRE
# a = ['a', 'a', 'a', 'a', 'b', 'b', 'c', 'c', 'c', 'd', 'e', 'e', 'e', 'e', 'e']

# letter_counts = Counter(features)
# df = pd.DataFrame.from_dict(letter_counts, orient='index')
# fig, ax = plt.subplots(figsize=(20, 10))
# df.plot(kind='bar', ax=ax)
# plt.title('Most Common Features Selected by MRMRE.')
# fig.savefig(f'features_{name}_{names[3]}.png' )
# #
# #
# print('Chosen in 80% of cases: \n', df[df[0] > 60])
# print('Chosen in 60% of cases: \n', df[df[0] > 60])
# print('Chosen in 40% of cases: \n', df[df[0] > 40])
# print('Chosen in 30% of cases: \n', df[df[0] > 30])
# print('Chosen in 20% of cases: \n', df[df[0] > 20])
