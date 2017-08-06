import pandas as pd
import argparse
"""
Takes in a csv, collects SMILES and molecule id information
and outputs a .smi file to pass to rfrag and indexing in rdkit

python make_smiles.py filename.csv 'column1,column2'

"""
parser = argparse.ArgumentParser()
parser.add_argument('filename', help = 'csv with smiles and molecule id')
parser.add_argument('columns', help = 'column headers separated by columns (ex: column1,column2)')
args = parser.parse_args()

filename = args.filename
columns = [col for col in args.columns.split(',')]

data = pd.read_csv(filename , header = 0)
smi = data[columns]
smi.to_csv(filename[:-3]+'smi', header = False, index = False)
