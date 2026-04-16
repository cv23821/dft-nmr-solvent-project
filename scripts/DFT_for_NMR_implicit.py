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
 
 
interval = 5
solvent = 'chloroform'
data = f'benchmark_{solvent}_optimized_xyz.txt'
# new_folder = 'run_orca/benchmark_chloroform_solvator_5_NMR_B3LYP'
 
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
    mol_name = mol_name_list[mol_index]
 
    # dft_method = 'WB97X'
    # basis_set = '6-311g(d,p)'
    # dft_method = 'B3LYP'
    # basis_set = 'cc-pVDZ'
    dft_method = 'PBE0'
    basis_set = 'def2-TZVP'
 
    out_folder = f'run_orca/benchmark_{solvent}_implicit_NMR_{dft_method}'
   
    out_folder = f'{out_folder}/{mol_name}'
    os.makedirs(out_folder, exist_ok=True)
   
    structure_file = f'run_orca/benchmark_{solvent}_implicit_optimized/{mol_name}/{mol_name}.xyz'
    with open(structure_file, 'r') as f:
        lines = f.read()
   
    rdmol = Chem.MolFromXYZBlock(lines)
    conf = rdmol.GetConformer()
    type_list = []
    xyz_list = []
   
    for atom in rdmol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        type_list.append(atom.GetSymbol())
        xyz_list.append([pos.x, pos.y, pos.z])

    rootline = f'{dft_method} {basis_set} D3 ALPB(CHCl3) TightSCF'
    write_block = write_orca_inp_block((type_list, xyz_list), rootline=rootline, control_block=None)
   
    write_block += '\n'
    write_block += f'%eprnmr\n'
    write_block += f' Nuclei = all C {{ shift, ssall}};\n'
    write_block += f' Nuclei = all N {{ shift, ssall}};\n'
    write_block += f' Nuclei = all F {{ shift, ssall}};\n'
    write_block += f' Nuclei = all Si {{ shift, ssall}};\n'
    write_block += f' Nuclei = all P {{ shift, ssall}};\n'
    write_block += f' Nuclei = all S {{ shift, ssall}};\n'
    write_block += f' Nuclei = all Cl {{ shift, ssall}};\n'
    write_block += f' Nuclei = all Br {{ shift, ssall}};\n'
    write_block += f' Nuclei = all H {{ shift, ssall}};\n'
    write_block += f' Nuclei = all O {{ shift, ssall}}\n'
    write_block += f' end\n'
   
    with contextlib.chdir(out_folder):
        inp_path = f'{mol_name}.inp'
        out_path = f'{mol_name}.out'
   
        with open(inp_path, 'w') as f:
            f.write(write_block)
        with open(out_path, 'w') as f:
            subprocess.run(['orca', inp_path], stdout=f, stderr=subprocess.STDOUT, check=False)
        

