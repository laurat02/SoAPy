"""
Example input script for generation of SLURM and Gaussian input files with geometries
obtained from user input.
"""

import SoAPy

import time
import os

cwd = os.getcwd()

# Molecule's geometry.
input_geometry = """ 
    1   O    0.000000    0.000000    0.000000
    2   C    0.000000    0.000000    1.408300
    3   C    1.448274    0.000000    1.938866
    4   O    2.106632   -1.228868    1.643301
    5   C    2.652189   -1.426647    0.371980
    6   O    1.645579   -1.579134   -0.610582
    7   C    3.553158   -0.269820   -0.100076
    8   O    4.726228   -0.317679    0.690477
    9   C    2.829686    1.083115    0.048487
    10  O    1.851727    1.225050   -0.949265
    11  C    2.254458    1.229486    1.478228
    12  O    3.317686    1.461156    2.374717
    13  H    0.422191    0.798390   -0.341631
    14  H   -0.536457   -0.899915    1.722221
    15  H   -0.526265    0.879326    1.813831
    16  H    1.402927    0.032196    3.038160
    17  H    3.252100   -2.350618    0.456761
    18  H    1.008432   -2.238039   -0.304286
    19  H    3.791432   -0.442172   -1.163566
    20  H    5.325112    0.398838    0.442517
    21  H    3.559337    1.895288   -0.092809
    22  H    1.414878    0.373313   -1.114565
    23  H    1.598609    2.109382    1.516559
    24  H    3.896690    0.679931    2.372625
    """


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
distance_threshold = ['0']
snapshots = ['1']
frequency = ['532']

t0 = time.time()

# What type of test is being researched, basis sets, distance thresholds, or number of snapshots?
dir_list, dir_parameters, relative_dir_list = SoAPy.functions.set_options(spectroscopy, shell_type, functional, basis_set, distance_threshold, snapshots, frequency)
t1 = time.time()

# Make the directory structure.
SoAPy.functions.make_directories(dir_list)

# Generate coordinate data from input geometry.
dir_data, num_solute_atoms, atom_types = SoAPy.get_geometry.coordinates_from_input(input_geometry, dir_list, dir_parameters)
t2 = time.time()

# Generate test specific input and SLURM files.
SoAPy.generate_files.make_cmpd_directory(dir_list, dir_parameters)
t3 = time.time()

# Modify the SLURM submission scripts.
SoAPy.generate_files.modify_SLURM(molecule_name, slurm_parameters, dir_list, dir_parameters)
t4 = time.time()

# Modify the input files.
SoAPy.generate_files.modify_input(molecule_name, gaussian_parameters, dir_list, dir_parameters, dir_data, num_solute_atoms, atom_types)
t5 = time.time()

# Generate batch submission script.
SoAPy.functions.generate_batch_submission_script(relative_dir_list, dir_parameters)
t6 = time.time()

print("Total Time: ", t6-t0)
