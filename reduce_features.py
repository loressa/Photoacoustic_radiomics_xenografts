import pandas
import csv
import numpy as np
import json
import pingouin as pg

# ----------------------------------------------------------------------------------------------------------------------
# Load value of features in csv with extra column for PatientName
file_features = 'All_features_corrected_final.csv'
features = pandas.read_csv(file_features)
selected_features = features.copy()

# ----------------------------------------------------------------------------------------------------------------------
# Load scores from make_table_classifiers.py
scores = json.load(open('dict_scores_final.json'))

# ----------------------------------------------------------------------------------------------------------------------
# Load list of discarded features
discarded = json.load(open('discard_kruskal_benjamini_0.25.json'))

# ----------------------------------------------------------------------------------------------------------------------
selected = 0
for key1 in features.keys():
    keep = True
    if (key1 in discarded):
        keep = False
    if ((key1 != 'PatientName') and (key1 != 'Model')):
        for key2 in features.keys():
            if ((key2 != 'PatientName') and (key2 != 'Model')):
                if (key1 != key2):
                    if (not key2 in discarded):
                        result_rm_corr = pg.rm_corr(data=features, x=key1, y=key2, subject='PatientName')
                        rm_corr = result_rm_corr.r.values[0]
                        if (rm_corr > 0.9): # Selected highly correlated with 90%
                            score1 = scores[key1]
                            score2 = scores[key2]
                            if (score2 > score1):
                                keep = False
                            if (score2 == score1):
                                print('There is a tie between %s and %s, choose one' %(key1, key2))
                                # This will choose the second one between the tie ones
                                if (key2 in selected_features.keys()):
                                    keep = False
    if (keep):
        selected = selected + 1
        print('Selected %s' %key1)
    else:
        selected_features.drop(columns=key1, inplace=True)

# Write output to csv
selected_features.to_csv('selected_final.csv')
print(selected)



