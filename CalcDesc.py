from rdkit import Chem
from rdkit.Chem import Descriptors, Descriptors3D, AllChem
import pandas as pd

"""
Calculate descriptors
"""


def calc_all(df_SMILES):
    '''
    get descriptors for all molecules loaded in from csv
    '''
    # Use nested dictionary to store calculated descriptors for each molecule
    # Create a dictionary
    d = {}
    def descriptor_calc(smiles, mol_name):
        '''
        Main function to calculate descriptors for molecules

        '''
        # read in molecules
        mol = (Chem.MolFromSmiles(smiles))

        # make parent dictionary
        d[mol_name]= {}
        # calculate descriptors and store as child dictionary
        for name, desc in Descriptors.descList:
             d[mol_name][name]= desc(mol)

        m2 = AllChem.AddHs(mol)
        AllChem.EmbedMolecule(m2)
        AllChem.MMFFOptimizeMolecule(m2)
        d[mol_name]['Asphericity'] = Descriptors3D.Asphericity(m2)
        d[mol_name]['PMI1'] = Descriptors3D.PMI1(m2)
        d[mol_name]['PMI2'] = Descriptors3D.PMI2(m2)
        d[mol_name]['PMI3'] = Descriptors3D.PMI3(m2)
        d[mol_name]['NPR1'] = Descriptors3D.NPR1(m2)
        d[mol_name]['NPR2'] = Descriptors3D.NPR2(m2)
        d[mol_name]['RadiusOfGyration'] = Descriptors3D.RadiusOfGyration(m2)
        d[mol_name]['InertialShapeFactor'] = Descriptors3D.InertialShapeFactor(m2)
        d[mol_name]['Eccentricity'] = Descriptors3D.Eccentricity(m2)
        d[mol_name]['SpherocityIndex'] = Descriptors3D.SpherocityIndex(m2)

    # -----------------------------
    for i in range(len(df_SMILES)):
        descriptor_calc(df_SMILES['SMILES'][i],
                                     df_SMILES['molecule'][i])

    return pd.DataFrame(d)
# --------------------------------------------------------------------------------
def add_IC50s(df_SMILES, df_desc):
    df_SMILES = df_SMILES.set_index('molecule')
    df_desc = df_desc.reset_index().T
    df_desc.columns = df_desc.iloc[0]
    df_desc = df_desc[1:]
    df_desc.index = df_desc.index.astype(int)
    return pd.merge(df_desc, df_SMILES, left_index = True, right_index = True)
# --------------------------------------------------------------------------------
def main(input_file, output_file, ic50s = False):
    df_SMILES = pd.read_csv(input_file)
    desc_df = calc_all(df_SMILES)

    if ic50s == True:
        desc_df = add_IC50s(df_SMILES, desc_df)
        desc_df.to_csv(output_file)
    else:
        desc_df.to_csv(output_file)
# --------------------------------------------------------------------------------
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 'Calculate descriptors')
    parser.add_argument('input_file', help = 'input file')
    parser.add_argument('output_file', help = 'output file name')
    parser.add_argument('-ic50s', action = 'store_true')
    args = parser.parse_args()

    if args.ic50s:
        main(args.input_file, args.output_file, ic50s = True)
    else:
        main(args.input_file, args.output_file)
