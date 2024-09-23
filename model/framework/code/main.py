# imports
import os
import csv
import sys
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]


root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(root)

from dilipred import DILIPRedictor

# read SMILES from .csv file, assuming one column with header
with open(input_file, "r") as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    smiles_list = [r[0] for r in reader]

R = []
dp = DILIPRedictor()
for s in smiles_list:
    r = dp.predict(s)
    r = r[["source", "value"]]
    r = r.set_index('source').T
    print(r)
    R += [r]


import pandas as pd

df = pd.concat(R)
df.to_csv(output_file, index=False)
