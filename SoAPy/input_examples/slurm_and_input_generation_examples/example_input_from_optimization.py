"""
Example input script for generation of SLURM and Gaussian input files with geometries
obtained from optimized structures. This script includes an optimization example.
"""

import SoAPy

import time
import os

cwd = os.getcwd()

# Location of optimized structures.
optimization_location = f'{cwd}'

# Molecule name.
molecule_name = 'H2O'

# Computational parameters.
slurm_parameters = {'time': '1-00:00',
                    'nodes': '1',
                    'cores': '24',
                    'queue': 'normal_q',
                    'allocation': 'your_allocation',
                    'email': 'your_email'}

gaussian_parameters = {'memory': '30GB',
                        'cores': '24',
                        'Opt': False,
                        'convergence': 'Tight',
                        'integrals': 'UltraFine',
                        'SCRF': False,
                        'solvent_model':'PCM',
                        'solvent': 'water',
                        'max_SCF_cycles': '512',
                        'Overlay': False,
                        'IOp': 'IOp(3/59=5)'} 


# Set testing parameters in lists.
spectroscopy = ['VCD', 'ROA']
shell_type = ['spherical']
functional = ['CAM-B3LYP']
basis_set = ['STO-3G']
distance_threshold = ['1.0']
snapshots = ['10']
frequency = ['532']

t0 = time.time()

# What type of test is being researched, basis sets, distance thresholds, or number of snapshots?
dir_list, dir_parameters, relative_dir_list = SoAPy.functions.set_options(spectroscopy, shell_type, functional, basis_set, distance_threshold, snapshots, frequency)
t1 = time.time()

# Make the directory structure.
SoAPy.functions.make_directories(dir_list)

# Generate coordinate data from input geometry.
dir_data, num_solute_atoms, atom_types = SoAPy.get_geometry.coordinates_from_optimization(relative_dir_list, dir_parameters, optimization_location)
t3 = time.time()

# Generate test specific input and SLURM files.
SoAPy.generate_files.make_cmpd_directory(dir_list, dir_parameters)
t4 = time.time()

# Modify the SLURM submission scripts.
SoAPy.generate_files.modify_SLURM(molecule_name, slurm_parameters, dir_list, dir_parameters)
t5 = time.time()

# Modify the input files.
SoAPy.generate_files.modify_input(molecule_name, gaussian_parameters, dir_list, dir_parameters, dir_data, num_solute_atoms, atom_types)
t6 = time.time()

# Generate batch submission script.
SoAPy.functions.generate_batch_submission_script(relative_dir_list, dir_parameters)
t7 = time.time()

print("Total Time: ", t7-t0)

