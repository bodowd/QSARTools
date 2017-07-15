from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np
from sklearn.externals import joblib
import requests

def fetch_WS(trgt):
    re = requests.get('https://www.ebi.ac.uk/chembl/api/data/target/{0}.json'.format(trgt))
    return (trgt, re.json()['pref_name'], re.json()['organism'])

def main(smiles, output, top = 10, model = '10uM'):
    if model == '1uM':
        morgan_nb = joblib.load('Data/models_23/1uM/mNB_1uM_all.pkl')
    else:
        morgan_nb = joblib.load('Data/models_23/10uM/mNB_10uM_all.pkl')

    classes = list(morgan_nb.targets)

    mol = Chem.MolFromSmiles(smiles)

    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits = 2048)
    res = np.zeros(len(fp), np.int32)
    DataStructs.ConvertToNumpyArray(fp, res)

    probas = list(morgan_nb.predict_proba(res.reshape(1,-1))[0])
    predictions = pd.DataFrame(list(zip(classes, probas)), columns=['id','probas'])

    top_pred = predictions.sort_values(by='probas', ascending = False).head(top)

    plist = []
    for i, e in enumerate(top_pred['id']):
        plist.append(fetch_WS(e))

    target_info = pd.DataFrame(plist, columns=['id', 'name', 'organism'])
    result = pd.merge(top_pred, target_info)
    result.to_csv(output)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('query', help = 'query molecule in SMILES')
    parser.add_argument('output', help = 'output filename')
    parser.add_argument('-model', help = '1uM or 10uM model', type = str)
    parser.add_argument('-top', help = 'top n hits', type = int)
    args = parser.parse_args()
    main(args.query, args.output, args.top, args.model)
