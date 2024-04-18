import SoAPy
import os
import pytest

def set_options_test():
    cwd = os.getcwd()

    # Molecule name
    molecule_name = 'TEST'

    # Set testing parameters in lists.
    spectroscopy = ['VCD', 'ROA', 'OptRot', 'Optimization']
    shell_type = ['molecular', 'spherical']
    functional = ['CAM-B3LYP', 'B3LYP']
    basis_set = ['6-31G(d)', 'cc-pVDZ', ['cc-pVDZ', '6-31G(d)'], 'Sadlej pVTZ', 'ORP', ]
    distance_threshold = ['2.5', '5']
    snapshots = ['10', '100']
    frequency = ['532', '489']

    # What type of test is being researched, basis sets, distance thresholds, or number of snapshots?
    dir_list, dir_parameters, relative_dir_list = SoAPy.functions.set_options(spectroscopy, shell_type, functional, basis_set, distance_threshold, snapshots, frequency)
    
    SoAPy.functions.make_directories(dir_list)
    
    for a in range(len(dir_list)):
        directory = os.path.exists(dir_list[a])
        assert directory == True

