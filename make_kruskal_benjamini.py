#!/usr/bin/env python

# Script written by Lorena Escudero (Department of Radiology, University of Cambridge)
# To compute Kruskal-Wallis p-values for a categorical variable (model) with Benjamini-Hochberg correction 

# Usage: Edit input file names and simply run as 'python3 make_kruskal_benjamini.py'

import scipy
from scipy import stats
import pandas

# OPT for plots
#from array import array
#from ROOT import TTree, TFile, TGraph, TGraphAsymmErrors, TCanvas, kRed, kBlue, TLegend, gStyle, TLine

# ---------------------------------------------------------------------------------------------------------------------

file_features = 'All_features_corrected_final.csv'
file_model = 'All_model.csv'

features = pandas.read_csv(file_features)
models = pandas.read_csv(file_model)

n_rows = len(features)
n_columns = len(features.keys())

# ---------------------------------------------------------------------------------------------------------------------   
# Calculation of Kruskal-Wallis
pval_list = []
feature_pvalue = {}
for key in features.keys():
    bin1_array, bin2_array = array('f'), array('f')
    for i in range(n_rows):
        if (models['Model'][i] == 'basal'):
            bin1_array.append(features[key][i])
        else:
            bin2_array.append(features[key][i])

    h_stat, p_value = stats.kruskal(bin1_array, bin2_array)
    pval_list.append(p_value)
    feature_pvalue[key] = p_value

# ---------------------------------------------------------------------------------------------------------------------   
# Benjamini-Hochberg correction
pval_list.sort()
Q = 0.25 # 25% false positives allowed
m = len(pval_list)
critical_values = []
for i in range(m):
    critical_values.append((i+1)*Q/m)

pval_critical = 100
for i in range(m):
    if (pval_list[i] >= critical_values[i]):
        if (pval_list[i] < pval_critical):
            pval_critical = pval_list[i]

# ---------------------------------------------------------------------------------------------------------------------   
# Compare p-values with critical one
counter_array, pval_array = array ('f'), array ('f')# for plots
counter_x = 0
for key in features.keys():
    counter_x = counter_x + 1
    p_value = feature_pvalue[key]
    if(p_value >= pval_critical):
        print('Discard %s' %key)
    counter_array.append(counter_x)
    pval_array.append(p_value)

# ---------------------------------------------------------------------------------------------------------------------   

# OPT: Make plots - using ROOT, so import if below is used
#graph = TGraph(n_columns, counter_array, pval_array)
#canvas = TCanvas("canvas", "canvas", 1000, 500)
#canvas.SetLeftMargin(0.15)
#canvas.SetRightMargin(0.15)
#canvas.SetTopMargin(0.15)
#canvas.SetBottomMargin(0.15)
#gStyle.SetOptStat(0)
#graph.SetTitle("")
#graph.SetMarkerStyle(20)
#graph.GetXaxis().SetTitle("Feature ID")
#graph.GetYaxis().SetTitle("Kruskal-Wallis p-value")
#graph.Draw('AP')
#line = TLine(0, pval_critical, 100, pval_critical)
#line.SetLineColor(kRed)
#line.Draw("same")

#name_figure = 'kruskal_pval_all_features_corrected_final.C'
#canvas.SaveAs(name_figure)
