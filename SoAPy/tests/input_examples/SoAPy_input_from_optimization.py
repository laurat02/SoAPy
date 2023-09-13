from SoAPy.functions import set_options
from SoAPy.functions import make_directories
from SoAPy.functions import generate_files
from SoAPy.get_geometry import coordinates_from_optimization
from SoAPy.functions import generate_batch_submission_script

import time
import os

cwd = os.getcwd()

# Molecule name.
molecule_name = 'R_fluorooxirane'

# Optimization location.
optimization_location = '/Users/bshumberger/Documents/Projects/basis_sets_project/SoAPy_tests/optimization/R_flourooxirane'

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
generate_files(molecule_name, dir_list, dir_parameters, dir_data, num_solute_atoms, atom_types)
t3 = time.time()

# Generate batch submission script.
generate_batch_submission_script(relative_dir_list, dir_parameters)

print("Total Time: ", t3-t0)

