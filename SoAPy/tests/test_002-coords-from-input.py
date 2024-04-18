import SoAPy
import os
import pytest

def test_coords_from_input():
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
    # Molecule name
    molecule_name = 'TEST'

    # Set testing parameters in lists.
    spectroscopy = ['VCD', 'ROA']
    shell_type = ['spherical']
    functional = ['CAM-B3LYP']
    basis_set = ['Sadlej pVTZ']
    distance_threshold = ['0']
    snapshots = ['1']
    frequency = ['532']

    # What type of test is being researched, basis sets, distance thresholds, or number of snapshots?
    dir_list, dir_parameters, relative_dir_list = SoAPy.functions.set_options(spectroscopy, shell_type, functional, basis_set, distance_threshold, snapshots, frequency)
    
    SoAPy.functions.make_directories(dir_list)
    
    dir_data, num_solute_atoms, atom_types = SoAPy.get_geometry.coordinates_from_input(input_geometry, dir_list, dir_parameters)
    
    with open(f"{cwd}/SoAPy/tests/test_002_coords_from_input/002ref.txt", "r") as file:
        for line in file:
            split_line = line.split("#")
        assert split_line[1] == str(dir_data)
        assert split_line[2] == str(num_solute_atoms)
        assert split_line[3] == str(atom_types)
        file.close()
    
