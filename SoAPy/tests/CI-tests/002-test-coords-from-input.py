import SoAPy
import pytest
import sys
import os

def test_coords_from_geom():
    cwd = os.getcwd()

    # Molecule name.
    molecule_name = 'R_fluorooxirane'

    # Molecule's geometry.
    input_geometry = """
        1   F   -1.4613   -0.0348    0.2295
        2   O    0.8156    0.7283    0.1265
        3   C   -0.2703    0.0081   -0.4705
        4   C    0.9160   -0.7016    0.1145
        5   H   -0.4252    0.0547   -1.5365
        6   H    0.7879   -1.2199    1.0554
        7   H    1.6306   -1.1467   -0.5647
        """

    # Set testing parameters in lists.
    spectroscopy = ['Optimization']
    shell_type = ['spherical']
    functional = ['CAM-B3LYP']
    basis_set = ['Sadlej pVTZ', 'ORP', 'LPol-ds', 'LPol-dl', 'LPol-fs', 'LPol-fl', 'aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ']
    distance_threshold = ['0']
    snapshots = ['1']
    frequency = []

    # What type of test is being researched, basis sets, distance thresholds, or number of snapshots?
    dir_list, dir_parameters, relative_dir_list = SoAPy.functions.set_options(spectroscopy, shell_type, functional, basis_set, distance_threshold, snapshots, frequency)

    # Make the directory structure.
    SoAPy.functions.make_directories(dir_list)

    # Generate coordinate data from input geometry.
    dir_data, num_solute_atoms, atom_types = SoAPy.coordinates_from_input(input_geometry, dir_list, dir_parameters)
    