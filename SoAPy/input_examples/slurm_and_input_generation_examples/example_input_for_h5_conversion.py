"""
Example input script for the generation of H5 file from a molecular dynamics trajectory.
"""

import SoAPy

import time
import os

cwd = os.getcwd()

# Location of MD trajectory.
trajectory_location = f'{cwd}/beta_D_glucose_example_trajectory.arc'

# Molecule name.
molecule_name = 'beta_D_glucose'

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
basis_set = ['Sadlej pVTZ']
distance_threshold = ['2.0']
snapshots = ['100']
frequency = ['532']

t0 = time.time()

# What type of test is being researched, basis sets, distance thresholds, or number of snapshots?
dir_list, dir_parameters, relative_dir_list = SoAPy.functions.set_options(spectroscopy, shell_type, functional, basis_set, distance_threshold, snapshots, frequency)
t1 = time.time()

# Make the directory structure.
SoAPy.functions.make_directories(dir_list)

# Make HDF5 file.
atoms, frames, hdf5_location = SoAPy.functions.convert_GROMACS(molecule_name, trajectory_location)

# Get atoms and frames from an already produced HDF5 file.
t2 = time.time()

print("Total Time: ", t7-t0)

