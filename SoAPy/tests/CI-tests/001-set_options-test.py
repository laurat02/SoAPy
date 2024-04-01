import SoAPy
import os

def test_set_options():
    cwd = os.getcwd()

    # Molecule name.
    molecule_name = 'Tester-001'

    # Set testing parameters in lists.
    spectroscopy = ['ROA', 'VCD']
    shell_type = ['spherical', ' molecular']
    functional = ['CAM-B3LYP', 'B3LYP']
    basis_set = ['Sadlej pVTZ', 'ORP', 'LPol-ds', 'LPol-dl', 'aug-cc-pVQZ', ['cc-pVDZ', 'STO-3G'], '6-31G(d)']
    distance_threshold = ['0', '1', '2']
    snapshots = ['10', '100']
    frequency = ['532']

    # What type of test is being researched, basis sets, distance thresholds, or number of snapshots?
    dir_list, dir_parameters, relative_dir_list = SoAPy.functions.set_options(spectroscopy, shell_type, functional, basis_set, distance_threshold, snapshots, frequency)

    # Make the directory structure.
    SoAPy.functions.make_directories(dir_list)
    
    for a in range(len(dir_list)):
        exists = os.path.exists(dir_list[a])
        assert exists == True
    