from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
# from rdkit.Chem import PandasTools
import pandas as pd
import json

df = pd.read_csv('data/JAK3_ChemBl.csv', header = 0)
mols = [Chem.MolFromSmiles(mol) for mol in df['SMILES']]
fps = [AllChem.GetMorganFingerprintAsBitVect(mol,2) for mol in mols]



nodes = []
edges = []
for i in range(len(fps[:100])):
    cmp_id = df['molecule'][i]
    smi = df['SMILES'][i]
    node = {'data':{'cmpd_id':cmp_id, 'smi':smi }}
    nodes.append(node)

for i in range(len(fps[:100])):
    for j in range(i):
        tc = DataStructs.TanimotoSimilarity(fps[i], fps[j])
        if tc >=0.6:
            source = df['molecule'][i]
            target = df['molecule'][j]
            edge = {'data':{'source':source, 'target':target}}
            edges.append(edge)

data = {'nodes':nodes, 'edges':edges}

f = open('mols.json', 'w')
f.write(json.dumps(data))
f.close()


