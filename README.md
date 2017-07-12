# QSARTools

### Helper tools for methods I use routinely 

#### CalcDesc.py will calculate descriptors using RDKitChem for molecules in csv format. It needs SMILES strings for the molecules
Command line tool for calculating descriptors:
python CalcDesc.py [-h] [-ic50s] input_file output_file
positional arguments:
  input_file  input file name
  output_file output file name
optional arguments:
-h, --help show this help message and exit
-ic50s  merge ic50s column to the newly created descriptor dataset





#### CV.py is a helper tool for cross validation using sklearn.
