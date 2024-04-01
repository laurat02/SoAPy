"""Provide the primary functions."""

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
import pandas as pd
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
                            #if num_functionals > 1:
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

            # Identify the snapshot separator.
            separator_identifier = content[0].split()
            separator_identifier = separator_identifier[1]

            #content.pop(0)

            # Iterate through the lines in the trajectory file.
            for i, line in enumerate(content):
                column = line.split()

                # Identifies the snapshot separator.
                #if column[1] != "Great":
                if column[1] != separator_identifier:
                    
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
                elif column[1] == separator_identifier or not line:
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



def get_atoms_frames_from_GROMACS(molecule_name, trajectory_location):
    """ 
    Obtains the atoms and frames assuming an HDF5 file has already been generated.
    """
    # Open trajectory.
    with open(f"{trajectory_location}","r") as traj:
        content = traj.readlines()
        lines = len(content)

        # Get number of atoms.
        atoms = int(content[0].split()[0])

        # Calculate number of frames.
        frames = int(lines / (atoms + 2))

    return atoms, frames



def generate_files(molecule_name, dir_list, dir_parameters, dir_data, num_solute_atoms, atom_types):
    """
    Creates a directory structure for a single test and modifies the input and SLURM files accordingly.
    """
    # Get current working directory.
    cwd = os.getcwd()

    # Custom basis sets of interest not currently set up in Gaussian but are in the basis set exchange.
    custom = ['ORP', 'Sadlej+', 'Sadlej pVTZ']

    # Custom basis sets of interest not currently in Gaussian or in the basis set exchange.
    custom1 = ['LPol-ds', 'LPol-dl', 'LPol-fs', 'LPol-fl']

    # Gets original input file and submission script locations.
    data_locations = os.path.join(os.path.dirname(getsourcefile(lambda:0)), "data/")
    input_file = os.path.join(data_locations, "input.dat")
    submission_file = os.path.join(data_locations, "G09_sub_SLURM.sh")
    
    # Setting up basis variable for pulling data from BSE.
    basis = 'X'

    # Change to test specific directory.
    for a in range(len(dir_list)):
        os.chdir(f"{dir_list[a]}")

        # Gets input file location if optimization is entered as the spectroscopy.
        if dir_parameters[a][0] == 'optimization':
            input_file = os.path.join(data_locations, "input_files/optimization/input.dat")

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
                if dir_parameters[a][0] != 'optimization':
                    content = content.replace("SPECTROSCOPY", dir_parameters[a][0])
                content = content.replace("FUNCTIONAL", dir_parameters[a][2])
                if len(dir_parameters[a][3]) != 2 and dir_parameters[a][3] not in custom and dir_parameters[a][3] not in custom1:
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

                # Pulling custom basis set data from the SoAPy basis set library.
                if dir_parameters[a][3] in custom1 and len(dir_parameters[a][3]) != 2:
                    current_basis = dir_parameters[a][3]
                    if current_basis != basis:
                        basis_file = os.path.join(data_locations, f"basis_sets/{dir_parameters[a][3]}.gbs")
                        basis_file = open(f"{basis_file}", "r")
                        lines = basis_file.readlines()
                        lines = np.asarray(lines)
                        for atom in atom_types:
                            atom_indices = np.where(lines == f"{atom}" + "  0\n" )
                            splt_indices = np.where(lines == "****\n")
                            for atom_index in range(len(atom_indices[0])):
                                for splt_index in range(len(splt_indices[0])):
                                    if splt_indices[0][splt_index] == atom_indices[0][atom_index] - 1:
                                        begin = atom_indices[0][atom_index]
                                        end = splt_indices[0][splt_index + 1]
                            count = begin
                            while count <= end:
                                file.write(lines[count])
                                count += 1

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

                    # Still need to include LPol basis sets in the mixed basis set section.

                
                #file.write('Surface=SAS')

                file.write("\n")
                file.write("\n")

            # Modifies the submission file with test specific data.
            with open(f"{dir_list[a]}/cmpd_{conformer_count}/G09_sub_SLURM.sh", "r+") as file:
                content = file.read()
                if len(dir_parameters[a][3]) != 2:
                    basis = dir_parameters[a][3].replace(" ", "")
                    content = content.replace("MOLECULE_NAME_CONFORMER_NUMBER", f"{molecule_name}_cmpd_{conformer_count}_{dir_parameters[a][0]}_{dir_parameters[a][1]}_shell_{dir_parameters[a][2]}_{basis}_{dir_parameters[a][4]}_{dir_parameters[a][5]}")
                elif len(dir_parameters[a][3]) == 2:
                    basis0 = dir_parameters[a][3][0].replace(" ", "")
                    basis1 = dir_parameters[a][3][1].replace(" ", "")
                    content = content.replace("MOLECULE_NAME_CONFORMER_NUMBER", f"{molecule_name}_cmpd_{conformer_count}_{dir_parameters[a][0]}_{dir_parameters[a][1]}_shell_{dir_parameters[a][2]}_{basis0}_{basis1}_{dir_parameters[a][4]}_{dir_parameters[a][5]}")
                file.seek(0)
                file.write(content)
                file.truncate()

            conformer_count += 1

    os.chdir(f"{cwd}")



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
    # Initialize directory level arrays.
    dir_frequencies = []
    dir_intensities = []
    dir_nbf = []

    # Change to test specific directory.
    for a in range(len(dir_list)):
        print("____________________________")
        print(f"Spectroscopy: {dir_parameters[a][0]} || Solvent Shell Type: {dir_parameters[a][1]} || Functional: {dir_parameters[a][2]} || Basis: {dir_parameters[a][3]} || Distance: {dir_parameters[a][4]} || Snapshots: {dir_parameters[a][5]}")
        with open(f"{cwd}/output_data.txt", "a") as file:
            file.write(f"Spectroscopy: {dir_parameters[a][0]} || Solvent Shell Type: {dir_parameters[a][1]} || Functional: {dir_parameters[a][2]} || Basis: {dir_parameters[a][3]} || Distance: {dir_parameters[a][4]} || Snapshots: {dir_parameters[a][5]}\n")

        test_frequencies = []
        test_intensities = []

        # Obtain data from each conformer in the test.
        conformer_count = 1 
        while conformer_count <= int(dir_parameters[a][5]):
            #print(f"Snapshot: {conformer_count}")
            frequency = []
            intensity = []
            num_imaginary_frequencies = None
            os.chdir(f"{dir_list[a]}/cmpd_{conformer_count}")

            # Data collection for ROA and VCD.
            if dir_parameters[a][0] != "OptRot":
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
                            if num_imaginary_frequencies == None:
                                num_imaginary_frequencies = 0
                            if split_line[0] == "NBasis=":
                                nbf = int(split_line[1])

                    # Correcting units according to the Gaussian output.
                    if dir_parameters[a][0] == "VCD":
                        intensity = [i * 10**-44 for i in intensity]
                    if dir_parameters[a][0] == "ROA":
                        intensity = [i * 10**4 for i in intensity]

                # Calculate total number of vibrations for nonlinear molecules.
                num_vibrations = 3 * natom - 6

                # Confirming that all vibrational frequencies have obtained from the Gaussian output file. This is required since some of the imaginary frequencies result in line splits that don't separate correctly.
                #if len(frequency) == len(intensity) and len(frequency) == num_vibrations:
                #    print("No discrepancies between frequencies and intensities obtained from output.")
                #else:
                #    print("Discrepancy between the frequncies, intensities, and number of vibrations.")

                # The possible discrepancy between line splits and the number of vibrational frequencies can be eliminated by removing the imaginary frequencies which we choose to do regardless of whether the descrepancy exists or not.
                real_frequencies = []
                real_intensities = []
                num_real_vibrations = num_vibrations - num_imaginary_frequencies
                j = num_vibrations - 1
                while j > num_imaginary_frequencies - 1:
                    real_frequencies.append(frequency[j])
                    real_intensities.append(intensity[j])
                    j -= 1

                print(f"GFE = {delta_G} \t Number of Atoms = {natom} \t Number of Basis Functions = {nbf} \t Vibrational Frequencies = {num_vibrations} \t Imaginary Frequencies = {num_imaginary_frequencies}")

                # Print frequencies and intensities to output file.
                with open(f"{cwd}/output_data.txt", "a") as file:
                    file.write(f"Snapshot = {conformer_count} \t GFE = {delta_G} \t Number of Atoms = {natom} \t Vibrational Frequencies = {num_vibrations} \t Imaginary Frequencies = {num_imaginary_frequencies}\n")
                    for i in range(num_real_vibrations):
                        file.write("{:.4f} \t {:e}\n".format(real_frequencies[i], real_intensities[i]))
                test_frequencies.extend(real_frequencies)
                test_intensities.extend(real_intensities)

            # Data collection for optical rotation.
            elif dir_parameters[a][0] == "OptRot":
                with open("output.log", "r") as file_out:
                    for line in file_out:
                        split_line = line.split()
                        if len(split_line) > 7:
                            if split_line[0] == "NAtoms=":
                                natom = int(split_line[1])
                            if split_line[0] == "Molar":
                                frequency.append(float(split_line[7]))
                                intensity.append(float(split_line[10]))
                            if split_line[0] == "NBasis=":
                                nbf = int(split_line[1])

                    # Correcting units on incident radiation back to nanometers.
                    frequency = [i / 10 for i in frequency]

                # Appending to test level arrays.
                test_frequencies.extend(frequency)
                test_intensities.extend(intensity)

                print(f"Number of Atoms = {natom} \t Number of Basis Functions = {nbf}")

                # Print frequencies and intensities to output file.
                with open(f"{cwd}/output_data.txt", "a") as file:
                    file.write(f"Snapshot = {conformer_count} \t Number of Atoms = {natom}\n")
                    for i in range(len(test_frequencies)):
                        file.write("{:.4f} \t {:e}\n".format(test_frequencies[i], test_intensities[i]))

            conformer_count += 1

        dir_frequencies.append(test_frequencies)
        dir_intensities.append(test_intensities)
        dir_nbf.append(nbf)

    return dir_frequencies, dir_intensities, dir_nbf



def generate_spectrum(fwhm, number_of_points, dir_list, dir_parameters, dir_frequencies, dir_intensities, min_frequency, max_frequency):
    """
    Generates the spectrum for each test. The spectral generation uses a simple averaging algorithm to weight the snapshots.
    """
    # Initialize directory level arrays for plotting.
    dir_frequency_axis = []
    dir_intensity_axis = []
    max_intensity = []

    for a in range(len(dir_list)):
        # Sorts the frequencies and intensities for a given test in ascending order of frequencies.
        spec_frequencies, spec_intensities = zip(*sorted(zip(dir_frequencies[a], dir_intensities[a])))
        
        # Define the interval at which points will be plotted for the x-coordinate.
        #delta = float((np.amax(np.array(dir_frequencies))-np.amin(np.array(dir_frequencies)))/number_of_points)
        delta = float((max_frequency - min_frequency)/number_of_points)

        # Compute the "spec_frequencies" array in electron volts (eV).
        spec_frequencies_eV = np.zeros_like(spec_frequencies)
        for i in range(len(spec_frequencies)):
            spec_frequencies_eV[i] = spec_frequencies[i]/8065.54429

        # Obtain the values associated with the x-coordinates in cm-1 and eV.        
        #frequency_axis = np.arange(np.amin(np.array(dir_frequencies)), np.amax(np.array(dir_frequencies)), delta)
        frequency_axis = np.arange(min_frequency, max_frequency, delta)
        frequency_axis_eV = frequency_axis/8065.54429

        # Initialize the array associated with the y-coordinates.
        intensity_axis = np.zeros_like(frequency_axis)

        # Normalize the intensity values based on the number of snapshots.
        normalized_spec_intensities = np.zeros_like(spec_intensities)
        for j in range(len(spec_intensities)):
            normalized_spec_intensities[j] = spec_intensities[j]*(1/int(dir_parameters[a][5]))

        # Fitting data to line shapes.
        # See equations 1.16 and 3.51 "Vibrational Optical Activity Principles and Applications" by Laurence Nafie for details.
        # Note the prefactor in equation 3.51 is excluded.
        if dir_parameters[a][0] == 'VCD':
            for b in range(len(frequency_axis)):
                for c in range(len(spec_frequencies)):
                    # Equation 8d in "ECD Cotton Effects Approximated by the Gaussian Curve and Other Methods" by Philip J. Stephens and Nobuyuki Harada.
                    #intensity_axis[b] += (1/((2.296*10**(-39))*np.sqrt(np.pi)*fwhm))*spec_frequencies_eV[c]*normalized_spec_intensities[c]*np.exp(-((frequency_axis_eV[b]-spec_frequencies_eV[c])/fwhm)**2)
                    intensity_axis[b] += spec_frequencies_eV[c]*normalized_spec_intensities[c]*(fwhm**2/(4*(frequency_axis_eV[b]-spec_frequencies_eV[c])**2+fwhm**2))
        if dir_parameters[a][0] == 'ROA':
            for b in range(len(frequency_axis)):
                for c in range(len(spec_frequencies)):
                    #intensity_axis[b] += normalized_spec_intensities[c]*(2/np.pi)*(fwhm/(4*(frequency_axis_eV[b]-spec_frequencies_eV[c])**2+fwhm**2))
                    intensity_axis[b] += spec_frequencies_eV[c]*normalized_spec_intensities[c]*(fwhm**2/(4*(frequency_axis_eV[b]-spec_frequencies_eV[c])**2+fwhm**2))

        # NEED TO CHECK UNITS FOR BOTH OF THESE EQUATIONS. THIS INCLUDES UNITS ON FWHM.
        # NEED TO ADD NORMALIZATION FUNCTION FOR THESE EQUATIONS.

        # Appending maximum intensity for normalization.
        max_intensity.append(max(np.absolute(intensity_axis)))

        # Appending frequencies and intensities to directory level arrays.
        dir_frequency_axis.append(frequency_axis)
        dir_intensity_axis.append(intensity_axis)

    return dir_frequency_axis, dir_intensity_axis, max_intensity
















