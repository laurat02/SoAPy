from SoAPy.functions import set_options
from SoAPy.functions import make_directories
from SoAPy.generate_files import make_cmpd_directory
from SoAPy.functions import convert_GROMACS
from SoAPy.functions import get_atoms_frames_from_GROMACS
from SoAPy.generate_files import modify_SLURM
from SoAPy.generate_files import modify_input
from SoAPy.get_geometry import coordinates_from_trajectory
from SoAPy.functions import generate_batch_submission_script

import time
import os

cwd = os.getcwd()

# Location of MD trajectory.
#trajectory_location = '/Users/bshumberger/Documents/SOlvation_Algorithm/SoAPy/SoAPy/tests/B_glucose_test_MD.arc'
#trajectory_location = '/Users/bshumberger/Documents/MD_trajectories/Water_6543.arc'
trajectory_location = '/Users/bshumberger/Documents/MD_trajectories/Water_2781.arc'

# Molecule name.
molecule_name = 'H2O'

# Computational parameters.
slurm_parameters = {'time': '6-00:00',
                    'nodes': '1',
                    'cores': '24',
                    'queue': 'normal_q',
                    'allocation': 'crawdad',
                    'email': 'bshumberger@vt.edu'}

gaussian_parameters = {'memory': '45GB',
                        'cores': '24',
                        'Opt': False,
                        'convergence': 'Tight',
                        'integrals': 'UltraFine',
                        'SCRF': False,
                        'solvent_model':'PCM',
                        'solvent': 'water',
                        'max_SCF_cycles': '512'} 


# Set testing parameters in lists.
spectroscopy = ['VCD', 'ROA']
shell_type = ['spherical']
functional = ['CAM-B3LYP']
basis_set = ['Sadlej pVTZ']
distance_threshold = ['2']
snapshots = ['100']
frequency = ['633nm']#, '589nm', '436nm', '355nm']

t0 = time.time()

# What type of test is being researched, basis sets, distance thresholds, or number of snapshots?
dir_list, dir_parameters, relative_dir_list = set_options(spectroscopy, shell_type, functional, basis_set, distance_threshold, snapshots, frequency)
t1 = time.time()

# Make the directory structure.
make_directories(dir_list)

# Make HDF5 file.
#atoms, frames, hdf5_location = convert_GROMACS(molecule_name, trajectory_location)

# Get atoms and frames from an already produced HDF5 file.
hdf5_location = '/Users/bshumberger/Documents/Projects/water_trajectories/H2O.h5'
atoms, frames = get_atoms_frames_from_GROMACS(molecule_name, trajectory_location)
t2 = time.time()

# Generate coordinate data from input geometry.
dir_data, num_solute_atoms, atom_types = coordinates_from_trajectory(atoms, frames, hdf5_location, dir_list, dir_parameters)
t3 = time.time()

# Generate test specific input and SLURM files.
make_cmpd_directory(dir_list, dir_parameters)
t4 = time.time()

# Modify the SLURM submission scripts.
modify_SLURM(molecule_name, slurm_parameters, dir_list, dir_parameters)
t5 = time.time()

# Modify the input files.
modify_input(molecule_name, gaussian_parameters, dir_list, dir_parameters, dir_data, num_solute_atoms, atom_types)
t6 = time.time()

# Generate batch submission script.
generate_batch_submission_script(relative_dir_list, dir_parameters)
t7 = time.time()

print("Total Time: ", t7-t0)

for a in range(0,len(dir_data[0])):
    print(dir_data[a])
