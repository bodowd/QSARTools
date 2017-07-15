# Predict drug targets for a query molecule based off ChemBl naive bayes classifier trained on ChemBl database

command line usage of ChemBl_pred.py

Example: 

python ChemBl_pred.py 'query smiles' 'output.csv' -model -top

-models  can be the 1uM or less stringent 10uM model

-top     int, returns that many predictions starting from highest probability

