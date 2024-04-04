from SoAPy.functions import set_options
from SoAPy.functions import convert_GROMACS
from SoAPy.functions import generate_files
from SoAPy.functions import generate_coordinates
from SoAPy.functions import generate_batch_submission_script
from SoAPy.functions import collect_data
import time
import os

cwd = os.getcwd()

# Location of MD trajectory.
trajectory_location = cwd + '/' + 'SoAPy/SoAPy/tests/B_glucose_test_MD.arc'

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

# Collect frequency and intensity data from Gaussian output files.
collect_data(cwd, dir_list, dir_parameters)
t2 = time.time()

print("Total Time: ", t2-t0)
