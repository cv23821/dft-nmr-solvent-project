import sys
import logging
import os
import argparse
import contextlib
import subprocess
import shutil
from rdkit import Chem


ems_path = '/user/home/cv23821/work'
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
solvent = 'water'
solvent_number = 15
data = f'lab_dataset_training_932_opt.txt'
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
    
    

    # dft_method = 'WB97X-D3'
    # basis_set = '6-311g(d,p)'
    # dft_method = 'B3LYP'
    # basis_set = 'cc-pVDZ'
    dft_method = 'PBE0'
    basis_set = 'def2-TZVP'
    out_folder = f'run_orca/lab_dataset_training_932_{solvent}_solvator_{solvent_number}_NMR_{dft_method}'
    out_folder = f'{out_folder}/{mol_name}'
    os.makedirs(out_folder, exist_ok=True)
    structure_file = f'run_orca/lab_dataset_training_932_water_solvator_15_optimized/{mol_name}/{mol_name}.xyz'
    with open(structure_file, 'r') as f:
        lines = f.read()
        num_atoms = int(lines.split('\n')[0].strip())
        num_atoms = num_atoms - solvent_number * 3    #change to num of atoms in solvent
    rdmol = Chem.MolFromXYZBlock(lines)
    conf = rdmol.GetConformer()
    type_list = []
    xyz_list = []
    for atom in rdmol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        type_list.append(atom.GetSymbol())
        xyz_list.append([pos.x, pos.y, pos.z])
    ###################
    H_idx_list = [str(i+1) for i, at in enumerate(type_list[:num_atoms]) if at == 'H']
    O_idx_list = [str(i+1) for i, at in enumerate(type_list[:num_atoms]) if at == 'O']
    S_idx_list = [str(i+1) for i, at in enumerate(type_list[:num_atoms]) if at == 'S']
    C_idx_list = [str(i+1) for i, at in enumerate(type_list[:num_atoms]) if at == 'C']
    Cl_idx_list = [str(i+1) for i, at in enumerate(type_list[:num_atoms]) if at == 'Cl']
    if len(H_idx_list) > 0:
        H_idx_string = ','.join(H_idx_list)
    if len(O_idx_list) > 0:
        O_idx_string = ','.join(O_idx_list)
    if len(S_idx_list) > 0:
        S_idx_string = ','.join(S_idx_list)
    if len(C_idx_list) > 0:
        C_idx_string = ','.join(C_idx_list)
    if len(Cl_idx_list) > 0:
        Cl_idx_string = ','.join(Cl_idx_list)
    rootline = f'{dft_method} {basis_set} D3 ALPB(WATER) TightSCF'
    write_block = write_orca_inp_block((type_list, xyz_list), rootline=rootline, control_block=None)
    write_block += '\n'
    write_block += f'%eprnmr\n'
    if len(C_idx_list) > 0:
        write_block += f' Nuclei = {C_idx_string} {{ shift, ssall}};\n'
    write_block += f' Nuclei = all N {{ shift, ssall}};\n'
    write_block += f' Nuclei = all F {{ shift, ssall}};\n'
    write_block += f' Nuclei = all Si {{ shift, ssall}};\n'
    write_block += f' Nuclei = all P {{ shift, ssall}};\n'
    if len(S_idx_list) > 0:
        write_block += f' Nuclei = {S_idx_string} {{ shift, ssall}};\n'
    if len(Cl_idx_list) > 0:
        write_block += f' Nuclei = {Cl_idx_string} {{ shift, ssall}};\n'
    write_block += f' Nuclei = all Br {{ shift, ssall}};\n'
    if len(H_idx_list) > 0:
        write_block += f' Nuclei = {H_idx_string} {{ shift, ssall}};\n'
    if len(O_idx_list) > 0:
        write_block += f' Nuclei = {O_idx_string} {{ shift, ssall}}\n'
    write_block += f' end\n'
    with contextlib.chdir(out_folder):
        inp_path = f'{mol_name}.inp'
        out_path = f'{mol_name}.out'
        with open(inp_path, 'w') as f:
            f.write(write_block)
        with open(out_path, 'w') as f:
            subprocess.run(['orca', inp_path], stdout=f, stderr=subprocess.STDOUT, check=False)

        # ==============================
        # CLEANUP IF ORCA SUCCESSFUL
        # ==============================
        if os.path.exists(out_path):
            with open(out_path, 'r') as f:
                content = f.read()

            if "ORCA TERMINATED NORMALLY" in content:
                logger.info(f"{mol_name} completed successfully — deleting scratch files.")

                deleted_files = 0
                for file in os.listdir('.'):
                    if file.endswith('.gbw') or file.startswith(f"{mol_name}.densities"):
                        os.remove(file)
                        deleted_files += 1

                logger.info(f"{deleted_files} scratch files deleted for {mol_name}.")
