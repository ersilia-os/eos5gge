#!/usr/bin/env python
# coding: utf-8

import os
import sys
from collections import Counter
import warnings

import numpy as np
import pandas as pd
from dimorphite_dl.dimorphite_dl import DimorphiteDL
from loguru import logger
from mordred import Calculator, descriptors
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, MolStandardize
from rdkit.Chem.MolStandardize import Standardizer, rdMolStandardize

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "code"))

from constants import (
    DESCS,
)

MODELPATH = os.path.join(root, "..", "..", "checkpoints")
FEATURESPATH = os.path.join(root, "..", "..", "checkpoints", "features")

DATA_PATH = os.path.join(root, "data", "drugbank_inchikeys.csv")

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

logger.remove()
logger.add(sys.stderr, level="CRITICAL")


def standardized_smiles(smiles):
    standardizer = Standardizer()

    # Read SMILES and convert it to RDKit mol object
    mol = Chem.MolFromSmiles(smiles)

    try:
        smiles_clean_counter = Counter()
        mol_dict = {}
        is_finalize = False

        for _ in range(5):

            # This solved phosphate oxidation in most cases but introduces a problem for some compounds: eg. geldanamycin where the stable strcutre is returned
            inchi_standardised = Chem.MolToInchi(mol)
            mol = Chem.MolFromInchi(inchi_standardised)

            # removeHs, disconnect metal atoms, normalize the molecule, reionize the molecule
            mol = rdMolStandardize.Cleanup(mol)
            # if many fragments, get the "parent" (the actual mol we are interested in)
            mol = rdMolStandardize.FragmentParent(mol)
            # try to neutralize molecule
            uncharger = (
                rdMolStandardize.Uncharger()
            )  # annoying, but necessary as no convenience method exists

            mol = uncharger.uncharge(mol)  # standardize molecules using MolVS and RDKit
            mol = standardizer.charge_parent(mol)
            mol = standardizer.isotope_parent(mol)
            mol = standardizer.stereo_parent(mol)

            # Normalize tautomers
            # Method 1
            normalizer = MolStandardize.tautomer.TautomerCanonicalizer()
            mol = normalizer.canonicalize(mol)

            # Method 2
            te = rdMolStandardize.TautomerEnumerator()  # idem
            mol = te.Canonicalize(mol)

            # Method 3
            mol = standardizer.tautomer_parent(mol)

            # Final Rules
            mol = standardizer.standardize(mol)
            mol_standardized = mol

            # convert mol object back to SMILES
            smiles_standardized = Chem.MolToSmiles(mol_standardized)

            if smiles == smiles_standardized:
                is_finalize = True
                break

            smiles_clean_counter[smiles_standardized] += 1
            if smiles_standardized not in mol_dict:
                mol_dict[smiles_standardized] = mol_standardized

            smiles = smiles_standardized
            mol = Chem.MolFromSmiles(smiles)

        if not is_finalize:
            # If the standardization process is not finalized, we choose the most common SMILES from the counter
            smiles_standardized = smiles_clean_counter.most_common()[0][0]
            # ... and the corresponding mol object
            # mol_standardized = mol_dict[smiles_standardized]

        return smiles_standardized

    except:

        return "Cannot_do"


def protonate_smiles(smiles):

    dimorphite = DimorphiteDL(min_ph=7.0, max_ph=7.0, pka_precision=0)
    protonated_smiles = dimorphite.protonate(smiles)

    # print("protonated_smiles")

    if len(protonated_smiles) > 0:
        protonated_smiles = protonated_smiles[0]

    return protonated_smiles


def smiles_to_inchikey(smiles):

    try:

        # Convert SMILES to a molecule object
        mol = Chem.MolFromSmiles(smiles)
        # Convert the molecule object to an InChI string
        inchi_string = Chem.MolToInchi(mol)
        # Convert the InChI string to an InChIKey
        inchi_key = Chem.inchi.InchiToInchiKey(inchi_string)

        return inchi_key

    except:

        return "Cannot_do"


def MorganFingerprint(s):
    x = Chem.MolFromSmiles(s)
    return AllChem.GetMorganFingerprintAsBitVect(x, 2, 2048)


def MACCSKeysFingerprint(s):
    x = Chem.MolFromSmiles(s)
    return AllChem.GetMACCSKeysFingerprint(x)


def get_num_charged_atoms_neg(mol):
    mol_h = Chem.AddHs(mol)
    Chem.rdPartialCharges.ComputeGasteigerCharges(mol_h)

    positive = 0
    negative = 0

    for atom in mol_h.GetAtoms():
        if float(atom.GetProp("_GasteigerCharge")) <= 0:
            negative += 1

    return negative


def get_num_charged_atoms_pos(mol):
    mol_h = Chem.AddHs(mol)
    Chem.rdPartialCharges.ComputeGasteigerCharges(mol_h)

    positive = 0
    # negative = 0

    for atom in mol_h.GetAtoms():
        if float(atom.GetProp("_GasteigerCharge")) >= 0:
            positive += 1
    return positive


def get_assembled_ring(mol):
    ring_info = mol.GetRingInfo()
    num_ring = ring_info.NumRings()
    ring_atoms = ring_info.AtomRings()
    num_assembled = 0

    for i in range(num_ring):
        for j in range(i + 1, num_ring):
            x = set(ring_atoms[i])
            y = set(ring_atoms[j])
            if not x.intersection(y):  # 2つの環が縮環でない場合に
                for x_id in x:
                    x_atom = mol.GetAtomWithIdx(x_id)
                    neighbors = [k.GetIdx() for k in x_atom.GetNeighbors()]
                    for x_n in neighbors:
                        if x_n in y:  # 環同士を繋ぐ結合があるか否か
                            num_assembled += 1

    return num_assembled


def get_num_stereocenters(mol):
    return AllChem.CalcNumAtomStereoCenters(
        mol
    ) + AllChem.CalcNumUnspecifiedAtomStereoCenters(mol)


def calc_descriptors(dataframe):
    mols = dataframe.smiles_r.apply(Chem.MolFromSmiles)
    # mols_fps=[AllChem.GetMorganFingerprintAsBitVect(x,2) for x in mols]
    descr = []
    for m in mols:
        descr.append(
            [
                Descriptors.TPSA(m),
                Descriptors.NumRotatableBonds(m),
                AllChem.CalcNumRings(m),
                Descriptors.NumAromaticRings(m),
                Descriptors.NumHAcceptors(m),
                Descriptors.NumHDonors(m),
                Descriptors.FractionCSP3(m),
                Descriptors.MolLogP(m),
                Descriptors.NHOHCount(m),
                Descriptors.NOCount(m),
                Descriptors.NumHeteroatoms(m),
                get_num_charged_atoms_pos(m),
                get_num_charged_atoms_neg(m),
                get_assembled_ring(m),
                get_num_stereocenters(m),
            ]
        )
    descr = np.asarray(descr)
    return descr


descs = DESCS


def calc_all_fp_desc(data):

    calc = Calculator(descriptors, ignore_3D=True)
    logger.debug(f"Calculated {len(calc.descriptors)} Descriptors")
    Ser_Mol = data["smiles_r"].apply(Chem.MolFromSmiles)
    # as pandas
    Mordred_table = calc.pandas(Ser_Mol)
    Mordred_table = Mordred_table.astype("float")
    # Mordred_table['smiles_r'] = model_tox_data['smiles_r']

    MACCSfingerprint_array = np.stack(data["smiles_r"].apply(MACCSKeysFingerprint))
    MACCS_collection = []
    for x in np.arange(MACCSfingerprint_array.shape[1]):
        x = "MACCS" + str(x)
        MACCS_collection.append(x)
    MACCSfingerprint_table = pd.DataFrame(
        MACCSfingerprint_array, columns=MACCS_collection
    )

    MorganFingerprint_array = np.stack(data["smiles_r"].apply(MorganFingerprint))
    Morgan_fingerprint_collection = []
    for x in np.arange(MorganFingerprint_array.shape[1]):
        x = "Mfp" + str(x)
        Morgan_fingerprint_collection.append(x)
    Morgan_fingerprint_table = pd.DataFrame(
        MorganFingerprint_array, columns=Morgan_fingerprint_collection
    )

    a = calc_descriptors(data)
    descdf = pd.DataFrame(a, columns=descs)
    descdf_approved = descdf.reset_index(drop=True)

    tox_model_data = pd.concat(
        [
            data,
            Morgan_fingerprint_table,
            MACCSfingerprint_table,
            descdf_approved,
            Mordred_table,
        ],
        axis=1,
    )

    # has_nans = tox_model_data.isna().any()
    # print(has_nans[has_nans == True])

    return tox_model_data

import csv
with open(DATA_PATH, "r") as f:
    reader = csv.reader(f)
    next(reader)
    smiles_list = []
    for r in reader:
        smiles_list += [r[0]]

data = {"smiles_r": smiles_list[:100]}
data = pd.DataFrame(data)

df = calc_all_fp_desc(data)
df = df[list(df.columns)[1:]]

all_na = df.columns[df.isna().all()]
df[all_na] = df[all_na].fillna(0)

from sklearn.impute import SimpleImputer

imputer = SimpleImputer(strategy="median")
imputer.fit(df)

import joblib

IMPUTER_FILE = os.path.join(root, "..", "..", "checkpoints", "imputer.joblib")
joblib.dump(imputer, IMPUTER_FILE)


print("Before", df.shape)
print("After", imputer.transform(df).shape)

print("Done!")