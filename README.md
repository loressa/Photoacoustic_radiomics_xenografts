Repository for the publication "Photoacoustic Image Radiomics in Patient-Derived Xenografts: A Study on Feature Sensitivity and Model Discrimination"

This repository contains:

*** DATA
- All_features_corrected_final.csv: spreadsheet with the values of the features, after correction for volume dependency has been applied.
- All_model.csv: a single column csv with the model corresponding to exactly the same order of columns as the file above.

*** SCRIPTS
- make_kruskal_benjamini.py: script that computes Kruskal-Wallis p-values, and the Benjamini-Hochberg correction.
- make_table_classifiers.py: script that creates a LaTeX table with the classification score per feature using different machine learning methods.
- reduce_features.py: script that runs the dimensionality reduction. 

Note: If input files of the scripts are not provided, it means they need to be created with some other script.

Pre-requisites/packages:
- csv
- json
- copy
- numpy
- pandas
- pingouin
- scipy
- sklearn
OPT: matplotlib and shap for plotting
