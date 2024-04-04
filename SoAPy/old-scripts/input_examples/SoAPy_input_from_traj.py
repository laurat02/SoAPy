from SoAPy.functions import set_options
from SoAPy.functions import make_directories
from SoAPy.functions import convert_GROMACS
from SoAPy.functions import generate_files
from SoAPy.get_geometry import coordinates_from_trajectory
from SoAPy.functions import generate_batch_submission_script

import time
import os

cwd = os.getcwd()

# Location of MD trajectory.
trajectory_location = '/Users/bshumberger/Documents/SOlvation_Algorithm/SoAPy/SoAPy/tests/B_glucose_test_MD.arc'

# Molecule name
molecule_name = 'B_glucose'

# Set testing parameters in lists.
spectroscopy = ['VCD']
shell_type = ['molecular']
functional = ['CAM-B3LYP']
basis_set = ['6-31G(d)', 'cc-pVDZ', ['cc-pVDZ', '6-31G(d)']]
distance_threshold = ['2.5']
snapshots = ['10']
frequency = []

t0 = time.time()

# What type of test is being researched, basis sets, distance thresholds, or number of snapshots?
dir_list, dir_parameters, relative_dir_list = set_options(spectroscopy, shell_type, functional, basis_set, distance_threshold, snapshots, frequency)
t1 = time.time()

# Make the directory structure.
make_directories(dir_list)

# Make HDF5 file
atoms, frames, hdf5_location = convert_GROMACS(molecule_name, trajectory_location)
t2 = time.time()

# Generate coordinate data for each required snapshot.
dir_data, num_solute_atoms, atom_types = coordinates_from_trajectory(atoms, frames, hdf5_location, dir_list, dir_parameters)
t3 = time.time()

# Generate test specific input and SLURM files.
generate_files(molecule_name, dir_list, dir_parameters, dir_data, num_solute_atoms, atom_types)
t4 = time.time()

# Generate batch submission script.
generate_batch_submission_script(relative_dir_list, dir_parameters)

print("Total Time: ", t4-t0)
