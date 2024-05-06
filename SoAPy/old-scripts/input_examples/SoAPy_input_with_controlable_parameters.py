from SoAPy.functions import set_options
from SoAPy.functions import make_directories
from SoAPy.generate_files import make_cmpd_directory
from SoAPy.generate_files import modify_SLURM
from SoAPy.generate_files import modify_input
from SoAPy.get_geometry import coordinates_from_optimization
from SoAPy.functions import generate_batch_submission_script

import time
import os

cwd = os.getcwd()

# Molecule name.
molecule_name = 'R_fluorooxirane'

# Optimization location.
optimization_location = '/Users/bshumberger/Documents/Projects/basis_sets_project/SoAPy_tests/optimization/R_flourooxirane'

# Computational parameters.
slurm_parameters = {'time': '6-00:00',
                    'nodes': '1',
                    'cores': '24',
                    'queue': 'normal_q',
                    'allocation': 'crawdad',
                    'email': 'bshumberger@vt.edu'}

gaussian_parameters = {'memory': '45GB',
                        'cores': '24',
                        'Opt': True,
                        'convergence': 'Tight',
                        'integrals': 'UltraFine',
                        'SCRF': False,
                        'solvent_model':'PCM',
                        'solvent': 'water',
                        'max_SCF_cycles': '512'} 


# Set testing parameters in lists.
spectroscopy = ['VCD']
shell_type = ['spherical']
functional = ['CAM-B3LYP']
basis_set = ['Sadlej pVTZ', 'ORP', 'LPol-ds', 'LPol-dl', 'LPol-fs', 'LPol-fl', 'aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ']
distance_threshold = ['0']
snapshots = ['1']
frequency = []

t0 = time.time()

# What type of test is being researched, basis sets, distance thresholds, or number of snapshots?
dir_list, dir_parameters, relative_dir_list = set_options(spectroscopy, shell_type, functional, basis_set, distance_threshold, snapshots, frequency)
t1 = time.time()

# Make the directory structure.
make_directories(dir_list)

# Generate coordinate data from input geometry.
dir_data, num_solute_atoms, atom_types = coordinates_from_optimization(relative_dir_list, dir_parameters, optimization_location)
t2 = time.time()

# Generate test specific input and SLURM files.
make_cmpd_directory(dir_list, dir_parameters)
t3 = time.time()

# Modify the SLURM submission scripts.
modify_SLURM(molecule_name, slurm_parameters, dir_list, dir_parameters)
t4 = time.time()

# Modify the input files.
modify_input(molecule_name, gaussian_parameters, dir_list, dir_parameters, dir_data, num_solute_atoms, atom_types)
t5 = time.time()

# Generate batch submission script.
generate_batch_submission_script(relative_dir_list, dir_parameters)
t6 = time.time()

print("Total Time: ", t6-t0)

