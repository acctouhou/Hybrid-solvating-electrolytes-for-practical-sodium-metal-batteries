
import numpy as np
from dataset.dataset import MoleculeDataset
import os
import shutil
import sys
import torch
import yaml
import numpy as np
from datetime import datetime

import torch.nn.functional as F
from torch.optim.lr_scheduler import CosineAnnealingLR

from utils.nt_xent import NTXentLoss

import pubchempy as pcp

import pandas as pd



#%%

data=pd.read_excel('data.xlsx')

#%%
code=data['code'].tolist()
name_list=data['name'].tolist()

    


#%%
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit.Chem import rdFreeSASA



log_1=[]
log_2=[]
log_3=[]
from tqdm import tqdm
import time
for i in tqdm(code):
    
    molecule = Chem.MolFromSmiles(i)
    
    
    molecule = Chem.AddHs(molecule)
    tpsa = Descriptors.TPSA(molecule)
    
    
    AllChem.EmbedMolecule(molecule,useRandomCoords=True)
    
    AllChem.UFFOptimizeMolecule(molecule)

    radii = rdFreeSASA.classifyAtoms(molecule)
    sasa = rdFreeSASA.CalcSASA(molecule, radii)


    log_1.append(tpsa)
    log_2.append(sasa)
    log_3.append(tpsa/sasa)
    


batch_size=128
temperature=0.1
use_cosine_similarity=True

device='cuda'
#%%
#from dataset.dataset_subgraph import MoleculeDataset,MoleculeDatasetWrapper
from dataset.dataset import MoleculeDataset,MoleculeDatasetWrapper

dataset=MoleculeDatasetWrapper(code).get_data_loaders()
from models.gcn_molclr import GCN
model = GCN(num_layer=5, emb_dim=300, feat_dim=512, drop_ratio=0, pool='mean').cuda()
model.load_state_dict(torch.load('model.pth',map_location='cuda:0'))
log=[]
with torch.no_grad():
    model.eval()
    valid_loss = 0.0
    counter = 0
    for xis in dataset:
        xis=xis.cuda()
        ris, zis = model(xis)
        log.append(ris.cpu().numpy())
x=np.vstack(log)
#%%
from joblib import dump, load
ox_clf=load('vox.joblib')
red_clf=load('vred.joblib')
dn_clf=load('dn.joblib')

vox=ox_clf.predict(x)
vre=red_clf.predict(x)

dn=dn_clf.predict(x)

#%%
import pickle


with open('result.pkl', 'rb') as file:
    esp_data = pickle.load(file)    
#%%
code=np.array(code)
min_esp=np.full_like(dn,np.nan)
max_esp=np.full_like(dn,np.nan)

for i in esp_data:
    
    temp1=str(i[0])
    temp2=i[1]
    if temp1 in code:
        sel=np.where(code==temp1)[0]
        min_esp[sel]=min(temp2)
        max_esp[sel]=max(temp2)
#%%
temp=np.loadtxt('homo.txt')        
homo=temp[:,0]
lumo=temp[:,1]




#%%
df = pd.DataFrame(np.array([name_list,min_esp,max_esp,
                            log_1,log_2,log_3,dn,vre,vox,homo,lumo]).T,
                  columns=['Name','ESP min(eV)','ESP max(eV)',
                           'Polar Surface Area(Å²)','Solvent Accessible Surface Area(Å²)',
                           'ratio','Donor(kcal/mol)', 'Vred V ( vs Li/Li+)', 'Vox V (vs Li/Li+)',
                           'HOMO(eV)','LUMO(eV)'])
df.to_excel('electrolyte properties_1006.xlsx')
