import sys
import logging
import os
import argparse
import contextlib
import subprocess

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
out_folder = 'run_orca/benchmark_water_solvator_5'
solvent_number = 15

atom_type_dict = {
    1: 'H',
    6: 'C',
    7: 'N',
    8: 'O',
    9: 'F',
    15: 'P',
    16: 'S',
    17: 'Cl',
    35: 'Br',
    53: 'I'}

# task_index = int(os.environ["SLURM_PROCID"])

parser = argparse.ArgumentParser()
parser.add_argument("bucket_index", type=int)
args = parser.parse_args()
bucket_index = args.bucket_index

running_index = bucket_index * interval

for i in range(interval):

    mol_index = running_index + i
    # if mol_index >= 932:
    #     break

    with open(f'datasets/{data}', 'r') as f:
        xyz_list = [line.strip() for line in f.readlines() if line.strip()]

    mol_file = xyz_list[mol_index]
    mol_name = mol_file.split('/')[-1].split('.')[0].split('_conf0')[0]

    emol = ems(mol_file, mol_id=mol_name)
    atom_types = [atom_type_dict[i] for i in emol.type]
    xyz = emol.xyz.tolist()

    rootline = 'GFN2-XTB ALPB(WATER)'
    control_block = f'%SOLVATOR\n NSOLV {solvent_number}\nEND'
    write_block = write_orca_inp_block((atom_types, xyz), rootline=rootline, control_block=control_block)

    orca_folder = f"{out_folder}/{mol_name}"
    os.makedirs(orca_folder, exist_ok=True)

    with contextlib.chdir(orca_folder):
        inp_path = f'{mol_name}.inp'
        out_path = f'{mol_name}.out'

        with open(inp_path, 'w') as f:
            f.write(write_block)
        with open(out_path, 'w') as f:
            subprocess.run(['orca', inp_path], stdout=f, stderr=subprocess.STDOUT, check=True)
