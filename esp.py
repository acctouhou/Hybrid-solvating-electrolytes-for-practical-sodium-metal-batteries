# Imports
import os
import argparse
import pickle
import numpy as np
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdMolDescriptors
from rdkit.Chem.rdForceFieldHelpers import UFFGetMoleculeForceField
from espsim.helpers import Renormalize, SimilarityMetric, psi4Charges, mlCharges, check_hs
import pandas as pd

# Define constants
methodPsi4 = 'B3LYP'
basisPsi4 = '6-311G(d,p)'
gridPsi4 = 32

# Argument parser setup (if needed)
parser = argparse.ArgumentParser()


# Load data
#loading the smiles code list

data = pd.read_excel('data.xlsx')['code'].tolist()

# Logging
log = []
log_fail = []

# Main processing loop
for i in tqdm(data):
    print(i)
    smiles = i
    try:
        # Convert SMILES to 3D structure
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.UFFOptimizeMolecule(mol)

        # Convert molecule to XYZ format
        xyz = Chem.rdmolfiles.MolToXYZBlock(mol)
        
        # Calculate charges using psi4
        charge = psi4Charges(xyz, basisPsi4, methodPsi4, gridPsi4, 'temp')
        print(charge)
        
        # Log results
        log.append([i, charge])
        
        # Rename output file for each molecule
        original_filename = '1_default_grid_esp.dat'
        new_filename = f'data/{smiles}.dat'
        os.rename(original_filename, new_filename)

        # Save log results to a pickle file
        with open('result.pkl', 'wb') as f:
            pickle.dump(log, f)

    except Exception as e:
        print(f"Error processing SMILES {smiles}: {e}")
        log_fail.append(smiles)
        np.save('temp', log_fail)

