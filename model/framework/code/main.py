# imports
import os
import csv
import sys
import pandas as pd
import numpy as np

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
sources = ['DILI', 'Diverse DILI C', 'BESP', 'Mitotox', 'Reactive Metabolite', 'Human hepatotoxicity', 'Animal hepatotoxicity A', 'Animal hepatotoxicity B', 'Preclinical hepatotoxicity', 'Diverse DILI A']
for s in smiles_list:
    try:
        r = dp.predict(s)
        r = r[["source", "value"]]
        r = r.set_index('source').T
        R += [r]
    except:
        r = pd.DataFrame([[np.nan]*len(sources)], columns=sources)
        R += [r]

# save results to .csv file
df = pd.concat(R)

with open(os.path.join(root, "..", "columns", "run_columns.csv"), "r") as f:
    reader = csv.reader(f)
    next(reader)
    columns = [r[0] for r in reader]

rename = {}
for i, c in enumerate(list(df.columns)):
    rename[c] = columns[i]

df = df.rename(columns=rename)
df.to_csv(output_file, index=False)
