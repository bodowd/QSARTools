import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('prefix', help = 'prefix to attach to file names')
args = parser.parse_args()

subprocess.check_call(["python ~/rdkit/Contrib/mmpa/rfrag.py < {}_smiles.smi > {}_frag.smi".format(args.prefix, args.prefix)], shell = True)
subprocess.check_call(["python ~/rdkit/Contrib/mmpa/indexing.py < {}_frag.smi > {}_mmp.csv".format(args.prefix, args.prefix)], shell = True)
