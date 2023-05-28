"""Provide the primary functions."""

import numpy as np
import h5py
import time 
import os
import shutil
import math
import requests
#import pandas as pd


from inspect import getsourcefile

def set_options(spectroscopy, shell_type, functional, basis_set, distance_threshold, snapshots, frequency=0):
    """ 
    Set the parameters and directory structure for testing.
    """
    # Get current working directory.
    cwd = os.getcwd()

    # Calculate the length of the lists of input parameters.
    num_spectroscopies = len(spectroscopy)
    num_shells = len(shell_type)
    num_functionals = len(functional)
    num_bases = len(basis_set)
    num_distances = len(distance_threshold)
    num_snapshots = len(snapshots)

    # Creates directory structure.
    relative_dir_list = []
    dir_list = []
    dir_parameters = []
    for a in range(num_spectroscopies):
        for b in range(num_shells):
            for c in range(num_functionals):
                for d in range(num_bases):
                    for e in range(num_distances):
                        for f in range(num_snapshots):
                            parameter_list = []
                            parameter_list.append(spectroscopy[a])
                            parameter_list.append(shell_type[b])
                            parameter_list.append(functional[c])
                            parameter_list.append(basis_set[d])
                            parameter_list.append(distance_threshold[e])
                            parameter_list.append(snapshots[f])
                            parameter_list.append(frequency)
                            if spectroscopy[a] == 'ROA':
                                parameter_list.append(frequency)
                            dir_name = spectroscopy[a]
                            if num_shells > 1:
                                dir_name = dir_name + '/' + shell_type[b]
                            if num_functionals > 1:
                                dir_name = dir_name + '/' + functional[c]
                            if num_bases > 1:
                                if len(basis_set[d]) == 2:
                                    dir_name = dir_name + '/' + basis_set[d][0] + '_' + basis_set[d][1]
                                else: 
                                    dir_name = dir_name + '/' + basis_set[d]
                            if num_distances > 1:
                                dir_name = dir_name + '/' + distance_threshold[e]
                            if num_snapshots > 1:
                                dir_name = dir_name + '/' + snapshots[f]
                            dir_list.append(cwd + '/' + dir_name)
                            relative_dir_list.append(dir_name)
                            dir_parameters.append(parameter_list)
    #for g in range(len(dir_list)):
        #os.makedirs(f"{dir_list[g]}")
    return dir_list, dir_parameters, relative_dir_list



def make_directories(dir_list):
    """
    Creates the directories for input generation.
    """
    for g in range(len(dir_list)):
        os.makedirs(f"{dir_list[g]}")



def convert_GROMACS(molecule_name, trajectory_location):
    """
    Converts a GROMACS MD trajectory file to an HDF5 file for python.
    """
    # Initialize arrays.
    atom_num = []
    atom_sym = []
    coord = []
    atom_type = []
    connect = []
    box_dims = np.zeros((3))
    
    frame_num = 1
    
    # Open trajectory file and write to HDF5 file.
    with h5py.File(f"{molecule_name}.h5", "w") as h5:
        with open(f"{trajectory_location}", "r") as traj:
            content = traj.readlines()
            lines = len(content)
            content.pop(0)
            
            # Iterate through the lines in the trajectory file.
            for i, line in enumerate(content):
                column = line.split()

                # Identifies the snapshot separator.
                if column[1] != "Great":

                    # Get box dimensions.
                    if column[-1] == "90.000000":
                        for a in range(3):
                            box_dims[a] = float(column[a])

                    # Get atom number, atomic symbol, coordinate data, atom type, and connectivity.
                    else:
                        atom_num.append(int(column[0]))

                        atom_sym.append(int(ord(str(column[1]))))

                        xyz_line = np.zeros((3))
                        for a in range(3):
                            xyz_line[a] = float(column[a+2])
                        coord.append(xyz_line)

                        atom_type.append(int(column[5]))

                        connect_line = np.zeros((4))
                        for a in range(len(column)-6):
                            connect_line[a] = int(column[6+a])
                        connect.append(connect_line)

                # Creates the datasets in the HDF5 file for the individual snapshot (group).
                elif column[1] == "Great" or not line:
                    grp = h5.create_group(f"frame{frame_num}")
                    grp.create_dataset("atom_num", data = np.array(atom_num))
                    grp.create_dataset('atom_sym', data = np.array(atom_sym))
                    grp.create_dataset("coord", data = np.array(coord))
                    grp.create_dataset("atom_type", data = np.array(atom_type))
                    grp.create_dataset("box_dims", data = np.array(box_dims))
                    grp.create_dataset("connect", data = np.array(connect))
                    frame_num += 1
                    atom_num = []
                    atom_sym = []
                    coord = []
                    atom_type = []
                    connect = []
                    box_dims = np.zeros((3))
            
            # Creates the last dataset in the HDF5 file for the individual snapshot (group).
            grp = h5.create_group(f"frame{frame_num}")
            grp.create_dataset("atom_num", data = np.array(atom_num))
            grp.create_dataset('atom_sym', data = np.array(atom_sym))
            grp.create_dataset("coord", data = np.array(coord))
            grp.create_dataset("atom_type", data = np.array(atom_type))
            grp.create_dataset("box_dims", data = np.array(box_dims))
            grp.create_dataset("connect", data = np.array(connect))   
     
    # Get number of atoms.
    with open(f"{trajectory_location}","r") as traj:
        atoms = int(traj.readline().split()[0])
 
    # Calculate number of frames.
    frames = int(lines / (atoms + 2))

    # Get location of HDF5 file.
    hdf5_location = os.path.join(os.getcwd(), f"{molecule_name}.h5")

    return atoms, frames, hdf5_location



def generate_files(molecule_name, dir_list, dir_parameters, dir_data, num_solute_atoms, atom_types):
    """
    Creates a directory structure for a single test and modifies the input and SLURM files accordingly.
    """
    # Get current working directory.
    cwd = os.getcwd()

    # Custom basis sets of interest not currently set up in Gaussian.
    custom = ['ORP', 'Sadlej+', 'Sadlej pVTZ']

    # Gets original input file and submission script locations.
    data_locations = os.path.join(os.path.dirname(getsourcefile(lambda:0)), "data/")
    input_file = os.path.join(data_locations, "input.dat")
    submission_file = os.path.join(data_locations, "G09_sub_SLURM.sh")
    
    # Setting up basis variable for pulling data from BSE.
    basis = 'X'

    # Change to test specific directory.
    for a in range(len(dir_list)):
        os.chdir(f"{dir_list[a]}")

        conformer_count = 1
        while conformer_count <= int(dir_parameters[a][5]):

            # Make directory for specific conformer.
            os.mkdir(f"cmpd_{conformer_count}")

            # Copies the original submission script into the conformer's directory.
            shutil.copy(f"{submission_file}", f"{dir_list[a]}/cmpd_{conformer_count}/")

            # Makes the submission script executable.
            os.chmod(f"{dir_list[a]}/cmpd_{conformer_count}/G09_sub_SLURM.sh", 0o755)

            # Copies the original input file into the conformer's directory.
            shutil.copy(f"{input_file}", f"{dir_list[a]}/cmpd_{conformer_count}/")

            # Modifies the input file with test specific data.
            with open(f"{dir_list[a]}/cmpd_{conformer_count}/input.dat", "r+") as file:
                content = file.read()
                content = content.replace("MOLECULE_NAME_CONFORMER_NUMBER", f"{molecule_name}_cmpd_{conformer_count}_{dir_parameters[a][0]}_{dir_parameters[a][1]}_shell")
                content = content.replace("SPECTROSCOPY", dir_parameters[a][0])
                content = content.replace("FUNCTIONAL", dir_parameters[a][2])
                if len(dir_parameters[a][3]) != 2 and dir_parameters[a][3] not in custom:
                    content = content.replace("BASIS_SET", dir_parameters[a][3])
                else:
                    content = content.replace("BASIS_SET", 'Gen')
                file.seek(0)
                file.write(content)
                file.truncate()
                
                # Writing atom and coordinate data from snapshot.
                for i in range(len(dir_data[a][conformer_count-1][0])):
                    file.write("{:<2} \t {:.6f} \t {:.6f} \t {:.6f}\n".format(dir_data[a][conformer_count-1][1][i], dir_data[a][conformer_count-1][2][i], dir_data[a][conformer_count-1][3][i], dir_data[a][conformer_count-1][4][i]))
                file.write("\n")

                # Writing frequency of incident radiation for ROA.
                if dir_parameters[a][0] == 'ROA':
                    for j in range(len(dir_parameters[a][6])):
                        file.write(f"{dir_parameters[a][6][j]}")
                        if len(dir_parameters[a][6]) != j + 1:
                            file.write(", ")
                        else:
                            file.write(" nm")
                    file.write("\n")
                    file.write("\n")

                # Pulling custom basis set data from the Basis Set Exchange API.
                if dir_parameters[a][3] in custom and len(dir_parameters[a][3]) != 2:
                    current_basis = dir_parameters[a][3]
                    if current_basis != basis:
                        main_bse_url = "http://basissetexchange.org"
                        base_url = os.environ.get('BSE_API_URL', main_bse_url)
                        headers = { 
                            'User-Agent': 'SOlvation Algorithm in PYthon (SoAPy)',
                            'From': 'bshumberger@vt.edu'
                        }
                        params = {'elements': atom_types}
                        custom_basis_data = requests.get(base_url + '/api/basis/' + dir_parameters[a][3] + '/format/gaussian94', params=params, headers=headers)
                        basis_data = custom_basis_data.text.split('!----------------------------------------------------------------------\n\n\n')
                        file.write(basis_data[1])
                        file.write("\n")
                        if custom_basis_data.status_code != 200:
                            raise RuntimeError("Could not obtain data from the BSE. Check the error information above.")
                        basis = current_basis
                    elif current_basis == basis:
                        file.write(basis_data[1])
                        file.write("\n")

                # Writing "Mixed" basis set data.
                if len(dir_parameters[a][3]) == 2:
                    file.write("1 - {} 0\n".format(num_solute_atoms))

                    if dir_parameters[a][3][0] in custom:
                        current_basis = dir_parameters[a][3][0]
                        if current_basis != basis:
                            main_bse_url = "http://basissetexchange.org"
                            base_url = os.environ.get('BSE_API_URL', main_bse_url)
                            headers = { 
                                'User-Agent': 'SOlvation Algorithm in PYthon (SoAPy)',
                                'From': 'bshumberger@vt.edu'
                            }   
                            params = {'elements': atom_types}
                            custom_basis_data = requests.get(base_url + '/api/basis/' + dir_parameters[a][3][0] + '/format/gaussian94', params=params, headers=headers)
                            basis_data = custom_basis_data.text.split('!----------------------------------------------------------------------\n\n\n')
                            file.write(basis_data[1])
                            if custom_basis_data.status_code != 200:
                                raise RuntimeError("Could not obtain data from the BSE. Check the error information above.")
                            basis = current_basis
                        elif current_basis == basis:
                            file.write(basis_data[1])
                    else:
                        file.write("{}\n".format(dir_parameters[a][3][0]))
                        file.write("****\n")

                    file.write("{} - {} 0\n".format(num_solute_atoms + 1, len(dir_data[a][conformer_count-1][0])))
                    if dir_parameters[a][3][1] in custom:
                        current_basis = dir_parameters[a][3][1]
                        if current_basis != basis:
                            main_bse_url = "http://basissetexchange.org"
                            base_url = os.environ.get('BSE_API_URL', main_bse_url)
                            headers = { 
                                'User-Agent': 'SOlvation Algorithm in PYthon (SoAPy)',
                                'From': 'bshumberger@vt.edu'
                            }   
                            params = {'elements': atom_types}
                            custom_basis_data = requests.get(base_url + '/api/basis/' + dir_parameters[a][3][1] + '/format/gaussian94', params=params, headers=headers)
                            basis_data = custom_basis_data.text.split('!----------------------------------------------------------------------\n\n\n')
                            file.write(basis_data[1])
                            file.write("\n")
                            if custom_basis_data.status_code != 200:
                                raise RuntimeError("Could not obtain data from the BSE. Check the error information above.")
                            basis = current_basis
                        elif current_basis == basis:
                            file.write(basis_data[1])
                            file.write("\n")
                    else:
                        file.write("{}\n".format(dir_parameters[a][3][1]))
                        file.write("****\n")
                        file.write("\n")

                file.write('Surface=SAS')

                file.write("\n")
                file.write("\n")

            # Modifies the submission file with test specific data.
            with open(f"{dir_list[a]}/cmpd_{conformer_count}/G09_sub_SLURM.sh", "r+") as file:
                content = file.read()
                if len(dir_parameters[a][3]) != 2:
                    content = content.replace("MOLECULE_NAME_CONFORMER_NUMBER", f"{molecule_name}_cmpd_{conformer_count}_{dir_parameters[a][0]}_{dir_parameters[a][1]}_shell_{dir_parameters[a][2]}_{dir_parameters[a][3]}_{dir_parameters[a][4]}_{dir_parameters[a][5]}")
                elif len(dir_parameters[a][3]) == 2:
                    content = content.replace("MOLECULE_NAME_CONFORMER_NUMBER", f"{molecule_name}_cmpd_{conformer_count}_{dir_parameters[a][0]}_{dir_parameters[a][1]}_shell_{dir_parameters[a][2]}_{dir_parameters[a][3][0]}_{dir_parameters[a][3][1]}_{dir_parameters[a][4]}_{dir_parameters[a][5]}")
                file.seek(0)
                file.write(content)
                file.truncate()

            conformer_count += 1

    os.chdir(f"{cwd}")

def generate_coordinates(atoms, frames, hdf5_location, dir_list, dir_parameters):
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
        for index in range(1, snapshots+1):
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

                atom_types = np.unique(input_sym)

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

    return dir_data, num_solute_atoms, atom_types



def generate_batch_submission_script(relative_dir_list, dir_parameters):
    """
    Creates a script to submit the jobs in batches.
    """
    # Get current working directory.
    cwd = os.getcwd()

    # Gets original batch submission script locations.
    data_locations = os.path.join(os.path.dirname(getsourcefile(lambda:0)), "data/")
    batch_submission_file = os.path.join(data_locations, "submit.py")

    # Copies the original submission script into the current directory.
    shutil.copy(f"{batch_submission_file}", f"{cwd}")

    # Inserts the list structures into the batch submision script.
    with open(f"{cwd}/submit.py", "r+") as file:
                content = file.read()
                content = content.replace("DIRECTORY_LIST", f"{relative_dir_list}")
                content = content.replace("DIRECTORY_PARAMETERS", f"{dir_parameters}")
                file.seek(0)
                file.write(content)
                file.truncate()



def collect_data(cwd, dir_list, dir_parameters):
    """
    Collects the data from the Gaussian output files for spectral generation.
    """
    # Change to test specific directory.
    for a in range(len(dir_list)):
        print("____________________________")
        print(f"Spectroscopy: {dir_parameters[a][0]} || Solvent Shell Type: {dir_parameters[a][1]} || Functional: {dir_parameters[a][2]} || Basis: {dir_parameters[a][3]} || Distance: {dir_parameters[a][4]} || Snapshots: {dir_parameters[a][5]}")
        with open(f"{cwd}/output_data.txt", "a") as file:
            file.write(f"Spectroscopy: {dir_parameters[a][0]} || Solvent Shell Type: {dir_parameters[a][1]} || Functional: {dir_parameters[a][2]} || Basis: {dir_parameters[a][3]} || Distance: {dir_parameters[a][4]} || Snapshots: {dir_parameters[a][5]}\n")

        # Obtain data from each conformer in the test.
        conformer_count = 1 
        while conformer_count <= int(dir_parameters[a][5]):
            print(f"Snapshot: {conformer_count}")
            frequency = []
            intensity = []
            os.chdir(f"{dir_list[a]}/cmpd_{conformer_count}")

            with open("output.log", "r") as file_out:
                for line in file_out:
                    split_line = line.split()
                    if len(split_line) > 7:
                        if split_line[0] == "Sum" and split_line[2] == "electronic" and split_line[6] == "Energies=":
                            delta_G = split_line[-1]
                    if len(split_line) > 2:
                        if split_line[2] == "imaginary":
                            num_imaginary_frequencies = int(split_line[1])
                        if split_line[0] == "NAtoms=":
                            natom = int(split_line[1])
                        if split_line[0] == "Frequencies":
                            frequency.extend(map(float, split_line[2:]))
                        if dir_parameters[a][0] == "VCD":
                            if split_line[0] == "Rot.":
                                intensity.extend(map(float, split_line[3:]))
                        elif dir_parameters[a][0] == "ROA":
                            if split_line[0] == "CID3":
                                intensity.extend(map(float, split_line[3:]))

            # Calculate total number of vibrations for nonlinear molecules.
            num_vibrations = 3 * natom - 6

            # Confirming that all vibrational frequencies have obtained from the Gaussian output file. This is required since some of the imaginary frequencies result in line splits that don't separate correctly.
            if len(frequency) == len(intensity) and len(frequency) == num_vibrations:
                print("No discrepancies between frequencies and intensities obtained from output.")
            else:
                print("Discrepancy between the frequncies, intensities, and number of vibrations.")
            # The possible discrepancy between line splits and the number of vibrational frequencies can be eliminated by removing the imaginary frequencies which we choose to do regardless of whether the descrepancy exists or not.
            real_frequencies = []
            real_intensities = []
            num_real_vibrations = num_vibrations - num_imaginary_frequencies
            j = num_vibrations - 1
            while j > num_imaginary_frequencies - 1:
                real_frequencies.append(frequency[j])
                real_intensities.append(intensity[j])
                j -= 1

            print(f"GFE = {delta_G} \t Number of Atoms = {natom} \t Vibrational Frequencies = {num_vibrations} \t Imaginary Frequencies = {num_imaginary_frequencies}")

            # Print frequencies and intensities to output file.
            with open(f"{cwd}/output_data.txt", "a") as file:
                file.write(f"Snapshot = {conformer_count} \t GFE = {delta_G} \t Number of Atoms = {natom} \t Vibrational Frequencies = {num_vibrations} \t Imaginary Frequencies = {num_imaginary_frequencies}\n")
                for i in range(num_real_vibrations):
                    file.write("{:.6f} \t \t {:.6f}\n".format(real_frequencies[i], real_intensities[i]))

            conformer_count += 1





















