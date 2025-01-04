# imports
import os
import csv
import sys
import pandas as pd

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
