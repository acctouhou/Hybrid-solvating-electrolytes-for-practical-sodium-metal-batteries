# Imports
import os
import argparse
import pickle
import numpy as np
import pandas as pd
import scipy.spatial
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdMolDescriptors
from rdkit.Chem.rdForceFieldHelpers import UFFGetMoleculeForceField
import psi4
from espsim.helpers import Renormalize, SimilarityMetric, psi4Charges, mlCharges, check_hs
import pandas as pd

# Define constants
methodPsi4 = 'B3LYP'
basisPsi4 = '6-311G(d,p)'
gridPsi4 = 32
psi4.set_num_threads(16) #enable multi-core

# Load data
#loading the smiles code list
data = pd.read_excel('data.xlsx')['code'].tolist()


# Logging
log = []
log_fail = []
count = 0

# Main processing loop
for i in tqdm(data):
    try:
        # Convert SMILES to 3D structure
        print(i)
        smiles = i
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.UFFOptimizeMolecule(mol)
        
        # Convert molecule to XYZ format and prepare psi4 geometry
        xyz = Chem.rdmolfiles.MolToXYZBlock(mol)
        psi4.geometry(xyz)
        
        # Set psi4 options
        psi4.set_options({
            'scf_type': 'df',
            'freeze_core': True,
            'basis': basisPsi4,
            'INTS_TOLERANCE': 1e-12,
            'e_convergence': 1e-8,
            'd_convergence': 1e-8,
        })
        
        # Calculate energy and extract HOMO-LUMO
        energy, wfn = psi4.energy(methodPsi4, return_wfn=True)
        homo = wfn.epsilon_a_subset("AO", "ALL").np[wfn.nalpha() - 1]
        lumo = wfn.epsilon_a_subset("AO", "ALL").np[wfn.nalpha()]
        
        # Append results to log
        print(homo, lumo)
        
        log.append([homo * 27.2114, lumo * 27.2114])
        #rescale to eV units
    except Exception as e:
        print(f"Error processing SMILES {i}: {e}")
        log.append([-100, -100])

# Save results to file
np.savetxt('homo.txt', np.array(log))