import sklearn
import pandas
import numpy as np
import copy

import sys
import json

from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier, GradientBoostingClassifier, AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn import svm

import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------------------------------------------------

# BDTs
max_depthValue = 3
n_estimatorsValue = 100
learning_rateValue = 1.0
random_stateValue = None

# SVM
kernelString = 'rbf'
kernelDegree = 2
gammaValue=0.05
coef0Value=1.0
cValue=1.0
tol=0.001
cache_size=1000

# ----------------------------------------------------------------------------------------------------------------------

file_features = 'All_features_corrected_final.csv'
file_model = 'All_model.csv'
dict_scores = {}

# ----------------------------------------------------------------------------------------------------------------------
# ALL FEATURES
features = pandas.read_csv(file_features)
print(features)
models = pandas.read_csv(file_model)

y = np.zeros(models.shape)
selected_model = 'luminal'
for i in range(models.shape[0]):
    if (models['Model'][i] == selected_model):
        y[i, 0] = 1

# ----------------------------------------------------------------------------------------------------------------------
# LATEX

original_stdout = sys.stdout
latex_file = 'table_classifiers_corrected_final.tex'
with open(latex_file, 'w') as latex:
    sys.stdout = latex
    print('\\documentclass[10pt]{article}')
    print('\\usepackage{xcolor}')
    print('\\usepackage{graphicx}')
    print('\\usepackage{authblk}')
    print('\\usepackage{chngpage}')
    print('\\begin{document}')
    print('\\begin{table}[ht]')
    print('\\footnotesize')
    print('\\hspace*{-3cm}')
    print('\\centering')
    print('\\begin{tabular}{|l|l|c|c|c|c|}')
    print('\\hline')
    print('ID & Feature & Random Forest & Gradient Boosting & Support Vector Machines & Average \\\\')
    #print('ID & Feature \\\\')
    print('\\hline')
    sys.stdout = original_stdout

# ----------------------------------------------------------------------------------------------------------------------
# FEATURES ONE BY ONE

n_rows = len(features)
n_columns = len(features.keys())

n_selected = 0
counter = 1
for key in features.keys():
    feature = np.zeros((n_rows, 1))
    for i in range(n_rows):
        feature[i] = features[key][i]

    rf = RandomForestClassifier(max_depth=max_depthValue, random_state=random_stateValue, n_estimators=n_estimatorsValue)
    gb = GradientBoostingClassifier(learning_rate=learning_rateValue, n_estimators=n_estimatorsValue, max_depth=max_depthValue)
    #bdt = AdaBoostClassifier(DecisionTreeClassifier(max_depth=max_depthValue),
                                  #n_estimators=n_estimatorsValue, learning_rate=learning_rateValue, random_state=random_stateValue)
    svm_model = svm.SVC(C=cValue, cache_size=cache_size, class_weight=None, coef0=coef0Value, decision_function_shape=None, degree=kernelDegree, gamma=gammaValue,
                       kernel=kernelString, max_iter=-1, probability=False, random_state=random_stateValue, shrinking=False, tol=tol, verbose=False)

    model_rf = rf.fit(feature, np.ravel(y, order='C'))
    model_gb = gb.fit(feature, np.ravel(y, order='C'))
    #model_bdt = bdt.fit(feature, np.ravel(y, order='C'))
    model_svm = svm_model.fit(feature, np.ravel(y, order='C'))

    score_rf = rf.score(feature, np.ravel(y, order='C'))
    score_gb = gb.score(feature, np.ravel(y, order='C'))
    #score_bdt = bdt.score(feature, np.ravel(y, order='C'))
    score_svm = svm_model.score(feature, np.ravel(y, order='C'))

    average_score = (score_rf + score_gb + score_svm) / 3.0
    dict_scores[key] = average_score

    with open(latex_file, 'a+') as latex:
        sys.stdout = latex
        row = str(counter) + ' & ' + key + ' & ' + ("%.3f" %score_rf) + ' & ' + ("%.3f" %score_gb) + ' & ' + ("%.3f" %score_svm) + ' & ' + ("%.3f" %average_score) + '\\\\'
        #row = str(counter) + ' & ' + key + '\\\\'
        print(row)
        print('\\hline')
        sys.stdout = original_stdout

    if (average_score > 0.7):
        print('Feature %s Score: %s' %(key, average_score))
        n_selected = n_selected + 1

    counter = counter + 1
print(n_selected)
json.dump(dict_scores,open('dict_scores_corrected_final.json','w'))

with open(latex_file, 'a+') as latex:
    sys.stdout = latex
    print('\\end{tabular}')
    print('\\end{table}')
    print('\\end{document}')
