"""Provides functions for obtaining geometries from different locations"""

import numpy as np
import h5py
import time 
import os
import shutil
import math
import requests
import matplotlib as mpl 
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from inspect import getsourcefile

def coordinates_from_trajectory(atoms, frames, hdf5_location, dir_list, dir_parameters):
    """
    Obtain snapshots from HDF5 file and modify coordinate data to center the solvent around the solute.
    """
    # Set solvent atom types.
    O = 349
    H = 350

    # Setting up data array.
    dir_data = []

    # Get snapshots for specified test.
    for a in range(len(dir_list)):
        test_data = []
        print(f'Spectroscopy: {dir_parameters[a][0]} || Solvent Shell Type: {dir_parameters[a][1]} || Functional: {dir_parameters[a][2]} || Basis: {dir_parameters[a][3]} || Distance: {dir_parameters[a][4]} || Snapshots: {dir_parameters[a][5]}')
        print(" ")
        snapshots = int(dir_parameters[a][5])
        flag = frames//int(snapshots)
        temp_index = flag

        # Set the file to be read.
        for index in range(1, snapshots + 1):
            snapshot_data = []
            print(f'Compound Number: {index}')
            with h5py.File(f"{hdf5_location}", "r") as h5:
                intermediate = f"frame{temp_index}/"
                atom_num = np.array(h5[intermediate]['atom_num'])
                atom_sym = np.array(h5[intermediate]['atom_sym'])
                coord = np.array(h5[intermediate]['coord'])
                atom_type = np.array(h5[intermediate]['atom_type'])
                box_dims = np.array(h5[intermediate]['box_dims'])
                connect = np.array(h5[intermediate]['connect'])

                # Initialize solute X, Y, and Z arrays.
                solute_num = []
                solute_sym = []
                solute_X = []
                solute_Y = []
                solute_Z = []

                # Initialize solvent X, Y, and Z arrays.
                solvent_num = []
                solvent_sym = []
                solvent_X = []
                solvent_Y = []
                solvent_Z = []

                # Appends solute and solvent arrays.
                for b in range(len(atom_num)):
                    if atom_type[b] != O and atom_type[b] != H:
                        solute_num.append(atom_num[b])
                        solute_sym.append(chr(atom_sym[b]))
                        solute_X.append(coord[b][0])
                        solute_Y.append(coord[b][1])
                        solute_Z.append(coord[b][2])
                    else:
                        solvent_num.append(atom_num[b])
                        solvent_sym.append(chr(atom_sym[b]))
                        solvent_X.append(coord[b][0])
                        solvent_Y.append(coord[b][1])
                        solvent_Z.append(coord[b][2])
                num_solute_atoms = len(solute_num)
                print(f'Number of Solute Atoms: {len(solute_num)}')
                print(f'Number of Solvent Atoms: {len(solvent_num)}')

                # If a solvent-only trajectory is provided, append the first solvent molecule to the solute arrays.
                if num_solute_atoms == 0:
                    solute_num.append(solvent_num[0])
                    solute_sym.append(solvent_sym[0])
                    solute_X.append(solvent_X[0])
                    solute_Y.append(solvent_Y[0])
                    solute_Z.append(solvent_Z[0])
                    
                    # Confirming presence of full solute molecule.
                    solute_complete = False
                    while not solute_complete:

                        # Initialize potential new atom arrays.
                        new_num = []
                        new_sym = []
                        new_X = []
                        new_Y = []
                        new_Z = []
                        for i in range(len(solute_num)):
                            atom_index = int(solute_num[i] - 1)
                            for j in range(len(connect[i])):

                                # Checking for duplicates.
                                if connect[atom_index][j] != 0:
                                    duplicate = False
                                    for k in range(len(solute_num)):
                                        if connect[atom_index][j] == solute_num[k]:
                                            duplicate = True
                                            break
                                        elif connect[atom_index][j] != solute_num[k]:
                                            for l in range(len(new_num)):
                                                if connect[atom_index][j] == new_num[l]:
                                                    duplicate = True
                                                    break
                                                elif connect[atom_index][j] != new_num[l]:
                                                    duplicate = False

                                    if duplicate == False:
                                        new_num.append(int(connect[atom_index][j]))
                                        solvent_index = int(connect[atom_index][j]) - 1

                                        new_sym.append(solvent_sym[solvent_index])
                                        new_X.append(solvent_X[solvent_index])
                                        new_Y.append(solvent_Y[solvent_index])
                                        new_Z.append(solvent_Z[solvent_index])

                        # Adding new atoms to the list which will be used in the Gaussian input file.
                        solute_num.extend(new_num)
                        solute_sym.extend(new_sym)
                        solute_X.extend(new_X)
                        solute_Y.extend(new_Y)
                        solute_Z.extend(new_Z)

                        # Checking to see if all solvent molecules are complete.
                        print("Number of New Solute Atoms: ", len(new_num))
                        if len(new_num) == 0:
                            solute_complete = True
                            print("Solute Complete")
                            break
                        else:
                            print("Solute Incomplete. Appending and obtaining connectivity of new atoms.")

                    # Remove previously appended solute atoms from solvent arrays.
                    sorted_solute = solute_num.copy()
                    sorted_solute.sort(reverse=True)
                    for m in sorted_solute:
                        del solvent_num[m-1]
                        del solvent_sym[m-1]
                        del solvent_X[m-1]
                        del solvent_Y[m-1]
                        del solvent_Z[m-1]

                # Compute average coordinate of solute molecule.
                avg_X = np.mean(solute_X)
                avg_Y = np.mean(solute_Y)
                avg_Z = np.mean(solute_Z)

                # Set conditions to center the solvent atoms around the solute.
                for c in range(len(solvent_num)):
                    if (solvent_X[c] - avg_X) > box_dims[0]/2:
                        solvent_X[c] -= box_dims[0]
                    elif (solvent_X[c] - avg_X) < -box_dims[0]/2:
                        solvent_X[c] += box_dims[0]
                    if (solvent_Y[c] - avg_Y) > box_dims[1]/2:
                        solvent_Y[c] -= box_dims[1]
                    elif (solvent_Y[c] - avg_Y) < -box_dims[1]/2:
                        solvent_Y[c] += box_dims[1]
                    if (solvent_Z[c] - avg_Z) > box_dims[2]/2:
                        solvent_Z[c] -= box_dims[2]
                    elif (solvent_Z[c] - avg_Z) < -box_dims[2]/2:
                        solvent_Z[c] += box_dims[2]

                # Computing spherical solvent shell arrays which are also used as intermediates for the molecular solvent shell arrays.
                # The spherical distance is defined as the distance of the farthest solute atom from the average coordinates of the solute plus the distance threshold.

                # Initializing spherical solvent shell arrays.
                sph_num = []
                sph_sym = []
                sph_X = []
                sph_Y = []
                sph_Z = []

                # Compute maximum distance of atoms in the solute from the average coordinate.
                max_dist = 0
                for d in range(len(solute_num)):
                    dist = math.sqrt((solute_X[d] - avg_X)**2 + (solute_Y[d] - avg_Y)**2 + (solute_Z[d] - avg_Z)**2)
                    if dist > max_dist:
                        max_dist = dist

                # Set the distance threshold for retaining solvent atoms.
                dist_thresh = max_dist + float(dir_parameters[a][4])

                # Compute the spherical solvent shell arrays.
                for e in range(len(solvent_num)):
                    sph_dist = math.sqrt((solvent_X[e] - avg_X)**2 + (solvent_Y[e] - avg_Y)**2 + (solvent_Z[e] - avg_Z)**2)
                    if sph_dist <= dist_thresh:
                        sph_num.append(solvent_num[e])
                        sph_sym.append(solvent_sym[e])
                        sph_X.append(solvent_X[e])
                        sph_Y.append(solvent_Y[e])
                        sph_Z.append(solvent_Z[e])

                print(f'Number of Spherical Solvent Atoms: {len(sph_num)}')

                # Computing the molecular solvent shell arrays.
                # The molecular solvent shell arrays are generated from the spherical solvent shell arrays to speed up the computations.

                if dir_parameters[a][1] == 'molecular':

                    # Initialize molecular solvent shell arrays.
                    mol_num = []
                    mol_sym = []
                    mol_X = []
                    mol_Y = []
                    mol_Z = []

                    # Determine distance of solvent atom from solute atoms and appends.
                    for f in range(len(solute_num)):
                        for g in range(len(sph_num)):
                            mol_dist = math.sqrt((solute_X[f] - sph_X[g])**2 + (solute_Y[f] - sph_Y[g])**2 + (solute_Z[f] - sph_Z[g])**2)
                            if mol_dist <= float(dir_parameters[a][4]):
                                no_duplicates = True
                                for h in range(len(mol_num)):
                                    if sph_num[g] == mol_num[h]:
                                        no_duplicates = False
                                        break
                                if no_duplicates:
                                    mol_num.append(sph_num[g])
                                    mol_sym.append(sph_sym[g])
                                    mol_X.append(sph_X[g])
                                    mol_Y.append(sph_Y[g])
                                    mol_Z.append(sph_Z[g])

                    print(f'Number of Atoms before Solvent Check: {len(mol_num)}')

            # Obtaining molecule data for the Gaussian input files.

            # Initialize arrays for Gaussian input files.
            input_num = []
            input_sym = []
            input_X = []
            input_Y = []
            input_Z = []

            # Appending solute data to Gaussian input arrays.
            input_num.extend(solute_num)
            input_sym.extend(solute_sym)
            input_X.extend(solute_X)
            input_Y.extend(solute_Y)
            input_Z.extend(solute_Z)

            # Appending solvent data to Gaussian input arrays based on solvent shell type.
            if dir_parameters[a][1] == 'spherical':
                input_num.extend(sph_num)
                input_sym.extend(sph_sym)
                input_X.extend(sph_X)
                input_Y.extend(sph_Y)
                input_Z.extend(sph_Z)
            elif dir_parameters[a][1] == 'molecular':
                input_num.extend(mol_num)
                input_sym.extend(mol_sym)
                input_X.extend(mol_X)
                input_Y.extend(mol_Y)
                input_Z.extend(mol_Z)

            # Confirming presence of full solvent molecule.
            solvent_complete = False
            while not solvent_complete:

                # Initialize potential new atom arrays.
                new_num = []
                new_sym = []
                new_X = []
                new_Y = []
                new_Z = []
                for i in range(len(input_num)):
                    atom_index = int(input_num[i] - 1)
                    for j in range(len(connect[i])):

                        # Checking for duplicates.
                        if connect[atom_index][j] != 0:
                            duplicate = False
                            for k in range(len(input_num)):
                                if connect[atom_index][j] == input_num[k]:
                                    duplicate = True
                                    break
                                elif connect[atom_index][j] != input_num[k]:
                                    for l in range(len(new_num)):
                                        if connect[atom_index][j] == new_num[l]:
                                            duplicate = True
                                            break
                                        elif connect[atom_index][j] != new_num[l]:
                                            duplicate = False

                            if duplicate == False:
                                new_num.append(int(connect[atom_index][j]))
                                solvent_index = int(connect[atom_index][j] - len(solute_num) - 1)

                                new_sym.append(solvent_sym[solvent_index])
                                new_X.append(solvent_X[solvent_index])
                                new_Y.append(solvent_Y[solvent_index])
                                new_Z.append(solvent_Z[solvent_index])

                # Adding new atoms to the list which will be used in the Gaussian input file.
                input_num.extend(new_num)
                input_sym.extend(new_sym)
                input_X.extend(new_X)
                input_Y.extend(new_Y)
                input_Z.extend(new_Z)

                # Checking to see if all solvent molecules are complete.
                print("Number of New Atoms: ", len(new_num))
                if len(new_num) == 0:
                    solvent_complete = True
                    print("Solvent Complete")
                    break
                else:
                    print("Solvent Incomplete. Appending and obtaining connectivity of new atoms.")

            all_types = np.unique(input_sym)

            ###
            solute_types = np.unique(solute_sym)
            solvent_types = np.unique(solvent_sym)
            atomic_symbol = {'H' : 1,
                             'He': 2,
                             'Li': 3,
                             'Be': 4,
                             'B' : 5,
                             'C' : 6,
                             'N' : 7,
                             'O' : 8,
                             'F' : 9}
            all_atom_types = sorted(all_types, key=lambda x: atomic_symbol[x])
            solute_atom_types = sorted(solute_types, key=lambda x: atomic_symbol[x])
            solvent_atom_types = sorted(solvent_types, key=lambda x: atomic_symbol[x])
            ###

            print("Total Number of Atoms: ", len(input_num))
            print(" ")

            # Combining snapshot specific data into one array.
            snapshot_data.append(input_num)
            snapshot_data.append(input_sym)
            snapshot_data.append(input_X)
            snapshot_data.append(input_Y)
            snapshot_data.append(input_Z)

            # Adding all data from each snapshot in a particular test to an array.
            test_data.append(snapshot_data)

            temp_index += flag

        # Adding all data to to be accessed later.
        dir_data.append(test_data)

    return dir_data, num_solute_atoms, all_atom_types, solute_atom_types, solvent_atom_types



def coordinates_from_input(input_geometry, dir_list, dir_parameters):
    """
    Generate coordinates from user provided input.
    This is currently set up to work with single molecules without explicit solvents.
    Additionally, the function can only read Cartesian coordinates as input.
    """
    # Read input geometry into array.
    geom = input_geometry.split()
    geometry = np.array(geom)

    # Getting number of atoms.
    num_solute_atoms = len(geometry) // 5

    # Separating coordinate data.
    input_num = []
    input_sym = []
    input_X = []
    input_Y = []
    input_Z = []

    for i in range(len(geometry)):
        if i % 5 == 0:
            input_num.append(int(geometry[i]))
        if i % 5 == 1:
            input_sym.append(str(geometry[i]))
        if i % 5 == 2:
            input_X.append(float(geometry[i]))
        if i % 5 == 3:
            input_Y.append(float(geometry[i]))
        if i % 5 == 4:
            input_Z.append(float(geometry[i]))

    # Combining snapshot specific data into one array.
    snapshot_data = []
    snapshot_data.append(input_num)
    snapshot_data.append(input_sym)
    snapshot_data.append(input_X)
    snapshot_data.append(input_Y)
    snapshot_data.append(input_Z)

    # Getting unique atoms.
    atom_types = np.unique(input_sym)

    # Setting up data array.
    dir_data = []

    # Setting up directory structure to access the coordinates.
    for a in range(len(dir_list)):
        test_data = []
        snapshots = int(dir_parameters[a][5])
        for index in range(1, snapshots+1):
            test_data.append(snapshot_data)
        dir_data.append(test_data)

    return dir_data, num_solute_atoms, atom_types



def coordinates_from_optimization(relative_dir_list, dir_parameters, optimization_location):
    """
    Reads geometries from Gaussian output files.
    This function assumes the same directory structure for the output data being obtained as the input data being produced.
    This function is primarily used to obtain optimized geometries for further testing.
    """
    # Dictionary of elements and their atomic numbers for transformation from Gaussian output.
    atomic_symbol = {1: 'H',
                    2: 'He',
                    3: 'Li',
                    4: 'Be',
                    5: 'B',
                    6: 'C',
                    7: 'N',
                    8: 'O',
                    9: 'F'}


    # Get current working directory.
    cwd = os.getcwd()

    # Setting up data array.
    dir_data = []

    # Go to test specific directory.
    for a in range(len(relative_dir_list)):
        test_data = []
        snapshots = int(dir_parameters[a][5])

        # Set the output file to be read from.
        for index in range(1, snapshots+1):
            output_file = f"{optimization_location}/{relative_dir_list[a]}/cmpd_{index}/output.log"
            if dir_parameters[a][0] == 'VCD':
                output_file = output_file.replace('VCD', 'optimization')
            if dir_parameters[a][0] == 'ROA':
                output_file = output_file.replace('ROA', 'optimization')

            # Get input orientation for the optimized geometry.
            with open(f"{output_file}", "r") as output:
                data = output.readlines()
                input_orientation_lines = []
                optimization_complete = False
            for line_num, line in enumerate(data):

                # Get number of atoms.
                if 'NAtoms' in line:
                    NAtom_line = line.split()
                    natoms = int(NAtom_line[1])
                    num_solute_atoms = natoms

                # Get all lines with input orientation. The last one is the optimized geometry.
                if 'Input orientation:' in line:
                    input_orientation_lines.append(line_num)

                # Confirm that a stationary point has been found.
                if 'Stationary point found.' in line:
                    optimization_complete = True
                if 'Number of steps exceeded' in line:
                    NStep_line = line.split()
                    nsteps = int(NStep_line[6])

            # Printing optimization status.
            print(f'Spectroscopy: {dir_parameters[a][0]} || Functional: {dir_parameters[a][2]} || Basis: {dir_parameters[a][3]} || Snapshots: {dir_parameters[a][5]}')            
            print('Optimization Complete: ', optimization_complete)
            if optimization_complete == False:
                print(f'Optimization terminated at step number {nsteps}.')

            # Get the beginning and ending lines of the geometry.
            begin = int(input_orientation_lines[-1] + 5)
            end = int(input_orientation_lines[-1] + 5 + natoms)

            # Separating coordinate data.
            input_num = []
            input_sym = []
            input_X = []
            input_Y = []
            input_Z = []

            for i in range(begin, end):
                geometry = data[i].split()
                input_num.append(int(geometry[0]))
                input_sym.append(str(atomic_symbol[int(geometry[1])]))
                input_X.append(float(geometry[3]))
                input_Y.append(float(geometry[4]))
                input_Z.append(float(geometry[5]))

            # Combining snapshot specific data into one array.
            snapshot_data = []
            snapshot_data.append(input_num)
            snapshot_data.append(input_sym)
            snapshot_data.append(input_X)
            snapshot_data.append(input_Y)
            snapshot_data.append(input_Z)

            # Getting unique atoms.
            atom_types = np.unique(input_sym)

            test_data.append(snapshot_data)

        dir_data.append(test_data)

    return dir_data, num_solute_atoms, atom_types


