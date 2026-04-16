import random
 
# path = 'datasets/lab_dataset_training.txt'
# output_file = 'datasets/lab_dataset_training_1000.txt'

 
# with open(path, 'r') as f:
#     content = f.readlines()
#     content = [line.strip() for line in content if line.strip()]
# random.shuffle(content)
# content = content[:1000]

 
# with open(output_file, 'w') as f:
#     for line in content:
#         f.write(line + '\n')


#path = '/home/b35bj/cv23821.b35bj/scratch/orca_test/run_orca/lab_dataset_training_1000_conformers'
#import os
# files = os.listdir(path)
# with open('datasets/lab_dataset_training_932_opt.txt', 'w') as f:
#     for file in files:
#         f.write(os.path.join(path, file) + '\n')

#this part is to get smile string of molecules

import os
import argparse
import sys

# ems_path = '/home/b35bj/cv23821.b35bj/scratch/'
# if ems_path not in sys.path:
#     sys.path.append(ems_path)

# from EMS.EMS import EMS as ems
# from EMS.modules.comp_chem.orca.orca_input import write_orca_inp_block

# path = 'orca_test/run_orca/lab_dataset_training_998_conformers_chcl3/XTB_ABELEZ_conf0.xyz'

# emol = ems(path)
# print(emol.mol_properties['SMILES'])

#code below is for naming

# import os

# path = '/home/b35bj/cv23821.b35bj/scratch/orca_test/run_orca/benchmark_water'

# mol_list = os.listdir(path)
# mol_list = [f'{path}/{i}' for i in mol_list]

# with open('/home/b35bj/cv23821.b35bj/scratch/orca_test/datasets/benchmark_water_optimized_xyz.txt', 'w') as f:
#     for mol in mol_list:
#         f.write(mol + '\n')

import sys
import logging
import os
import argparse
import contextlib
import subprocess
import shutil

from rdkit import Chem
from rdkit.Chem import Draw


smiles = 'O=C(NCCCCC(N)C(=O)O)C'
rdmol = Chem.MolFromSmiles(smiles, sanitize=True)
 
for atom in rdmol.GetAtoms():
    atom.SetProp("atomLabel", str(atom.GetIdx()))
 
img = Draw.MolToImage(rdmol, size=(400, 400))
img.save("molecule_with_indices.jpg")
