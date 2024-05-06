from SoAPy.functions import set_options
from SoAPy.functions import collect_data
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
spectroscopy = ['OptRot', 'ROA']
shell_type = ['spherical']
functional = ['CAM-B3LYP']
basis_set = ['6-31G(d)']
distance_threshold = ['2']
snapshots = ['1']
frequency = ['633nm']#, '589nm', '436nm', '355nm']

t0 = time.time()

# What type of test is being researched, basis sets, distance thresholds, or number of snapshots?
dir_list, dir_parameters, relative_dir_list = set_options(spectroscopy, shell_type, functional, basis_set, distance_threshold, snapshots, frequency)
t1 = time.time()

# Collect frequency and intensity data from Gaussian output files.
dir_frequencies, dir_intensities, dir_nbf = collect_data(cwd, dir_list, dir_parameters)
t2 = time.time()

print("Total Time: ", t2-t0)

