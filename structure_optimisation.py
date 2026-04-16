import sys
import logging
import os
import argparse
import contextlib
import subprocess
import shutil

from rdkit import Chem

ems_path = '/home/b35bj/cv23821.b35bj/scratch'
if ems_path not in sys.path:
    sys.path.append(ems_path)

from EMS.EMS import EMS as ems
from EMS.modules.comp_chem.orca.orca_input import write_orca_inp_block

########### Set up the logger system ###########
logger = logging.getLogger(__name__)
stdout = logging.StreamHandler(stream = sys.stdout)
formatter = logging.Formatter("%(name)s: %(asctime)s | %(levelname)s | %(filename)s:%(lineno)s >>> %(message)s")
stdout.setFormatter(formatter)
logger.addHandler(stdout)
logger.setLevel(logging.INFO)
########### Set up the logger system ###########


interval = 1
data = 'benchmark_water_optimized_xyz.txt'
# out_folder = 'run_orca/lab_dataset_training_932_opt_solvator_20dmso'
new_folder = 'run_orca/benchmark_water_solvator_15_optimized'



with open(f'datasets/{data}', 'r') as f:
    mol_name_list = [line.strip() for line in f.readlines() if line.strip()]
    mol_name_list = [name.split('/')[-1].split('.')[0].split('_conf0')[0] for name in mol_name_list]


parser = argparse.ArgumentParser()
parser.add_argument("bucket_index", type=int)
args = parser.parse_args()
bucket_index = args.bucket_index

run_index = bucket_index * interval


for i in range(interval):
    
    mol_index = run_index + i


    mol2_name = mol_name_list[mol_index]

    path = '/home/b35bj/cv23821.b35bj/scratch/orca_test/run_orca/benchmark_water_solvator_15'
    full_path = f'{path}/{mol2_name}/{mol2_name}.solvator.xyz'
    
    new_folder = f'{new_folder}/{mol2_name}'
    
    os.makedirs(new_folder, exist_ok=True)
    new_file = f'{new_folder}/{mol2_name}.solvator.xyz'
    shutil.copy(full_path, new_file)
    
    with open(new_file, 'r') as f:
        lines = f.read()
        print(lines)
    
    rdmol = Chem.MolFromXYZBlock(lines)
    
    conf = rdmol.GetConformer()
    type_list = []
    xyz_list = []
    
    for atom in rdmol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
    
        type_list.append(atom.GetSymbol())
        xyz_list.append([pos.x, pos.y, pos.z])
    
    rootline = 'GFN2-xTB ALPB(H2O) OPT TightSCF'
    control_block = None
    write_block = write_orca_inp_block((type_list, xyz_list), rootline=rootline, control_block=control_block)
    
    with contextlib.chdir(new_folder):
        inp_path = f'{mol2_name}.inp'
        out_path = f'{mol2_name}.out'
    
        with open(inp_path, 'w') as f:
            f.write(write_block)
        with open(out_path, 'w') as f:
            subprocess.run(['orca', inp_path], stdout=f, stderr=subprocess.STDOUT)