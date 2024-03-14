"""Provides functions associated with creating and modifying SLURM scripts and input files."""

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



def make_cmpd_directory(dir_list, dir_parameters):
    """
    Makes the compound (snapshot) level directories and copy SLURM script and input file.
    """
    # Get current working directory.
    cwd = os.getcwd()

    # Gets original input file and submission script locations.
    data_locations = os.path.join(os.path.dirname(getsourcefile(lambda:0)), "data/submission_and_input_files/")
    input_file = os.path.join(data_locations, "input.dat")
    submission_file = os.path.join(data_locations, "G09_sub_SLURM.sh")

    # Change to test specific directory.
    for a in range(len(dir_list)):
        os.chdir(f"{dir_list[a]}")
        snapshots = int(dir_parameters[a][5])

        for conformer_count in range(1, snapshots+1): 

            # Make directory for specific conformer.
            os.mkdir(f"cmpd_{conformer_count}")

            # Copies the original submission script into the conformer's directory.
            shutil.copy(f"{submission_file}", f"{dir_list[a]}/cmpd_{conformer_count}/")

            # Makes the submission script executable.
            os.chmod(f"{dir_list[a]}/cmpd_{conformer_count}/G09_sub_SLURM.sh", 0o755)

            # Copies the original input file into the conformer's directory.
            shutil.copy(f"{input_file}", f"{dir_list[a]}/cmpd_{conformer_count}/")

    os.chdir(f"{cwd}")



def modify_SLURM(molecule_name, slurm_parameters, dir_list, dir_parameters):
    """
    Modifies the specified SLURM file.
    """
    # Sets up loop structure.
    for a in range(len(dir_list)):
        snapshots = int(dir_parameters[a][5])
        for conformer_count in range(1, snapshots+1):

            # Opens the SLURM file being modified.
            with open(f"{dir_list[a]}/cmpd_{conformer_count}/G09_sub_SLURM.sh", "r+") as file:
                content = file.read()

                # Modifies the spacing in the basis to accomadate the submission file format.
                if len(dir_parameters[a][3]) != 2:
                    basis = dir_parameters[a][3].replace(" ", "") 
                    content = content.replace("MOLECULE_NAME_CONFORMER_NUMBER", f"{molecule_name}_cmpd_{conformer_count}_{dir_parameters[a][0]}_{dir_parameters[a][1]}_shell_{dir_parameters[a][2]}_{basis}_{dir_parameters[a][4]}_{dir_parameters[a][5]}")
                elif len(dir_parameters[a][3]) == 2:
                    basis0 = dir_parameters[a][3][0].replace(" ", "") 
                    basis1 = dir_parameters[a][3][1].replace(" ", "") 
                    content = content.replace("MOLECULE_NAME_CONFORMER_NUMBER", f"{molecule_name}_cmpd_{conformer_count}_{dir_parameters[a][0]}_{dir_parameters[a][1]}_shell_{dir_parameters[a][2]}_{basis0}_{basis1}_{dir_parameters[a][4]}_{dir_parameters[a][5]}")

                # Replaces keywords with SLURM parameters.
                content = content.replace("time", slurm_parameters["time"])
                content = content.replace("nodes", slurm_parameters["nodes"])
                content = content.replace("cores", slurm_parameters["cores"])
                content = content.replace("queue", slurm_parameters["queue"])
                content = content.replace("allocation", slurm_parameters["allocation"])
                content = content.replace("email", slurm_parameters["email"])
                file.seek(0)
                file.write(content)
                file.truncate()



def modify_input(molecule_name, gaussian_parameters, dir_list, dir_parameters, dir_data, num_solute_atoms, all_atom_types, solute_atom_types=0, solvent_atom_types=0):
    """
    Modifies the specified input files.
    """
    # Gets original input file and submission script locations.
    data_locations = os.path.join(os.path.dirname(getsourcefile(lambda:0)), "data/")

    # Custom basis sets of interest not currently set up in Gaussian but are in the basis set exchange.
    custom = ['ORP', 'Sadlej+', 'Sadlej pVTZ']

    # Custom basis sets of interest not currently in Gaussian or in the basis set exchange.
    custom1 = ['LPol-ds', 'LPol-dl', 'LPol-fs', 'LPol-fl', 'augD-3-21G', 'augT3-3-21G', 'R-ORP', 'rDPS']

    # Setting up basis variable for pulling data from BSE.
    basis = 'X'

    # Sets up loop structure.
    for a in range(len(dir_list)):
        snapshots = int(dir_parameters[a][5])
        for conformer_count in range(1, snapshots+1):

            # Opens the input file being modified.
            with open(f"{dir_list[a]}/cmpd_{conformer_count}/input.dat", "r+") as file:
                content = file.read()

                # Replaces keywords with Gaussian parameters.
                content = content.replace("MOLECULE_NAME_CONFORMER_NUMBER", f"{molecule_name}_cmpd_{conformer_count}_{dir_parameters[a][1]}_shell_{dir_parameters[a][3]}")
                content = content.replace("checkpoint_file", f"{molecule_name}_cmpd_{conformer_count}")
                content = content.replace("memory", gaussian_parameters["memory"])
                content = content.replace("cores", gaussian_parameters["cores"])
                content = content.replace("FUNCTIONAL", dir_parameters[a][2])

                if gaussian_parameters["Opt"]:
                    content = content.replace("convergence", gaussian_parameters["convergence"])
                else:
                    content = content.replace("Opt=convergence ", "")
                content = content.replace("integrals", gaussian_parameters["integrals"])
                if gaussian_parameters["SCRF"]:
                    content = content.replace("solvent_model", gaussian_parameters["solvent_model"])
                    content = content.replace("solvent", gaussian_parameters["solvent"])
                else:
                    content = content.replace("SCRF=(solvent_model, Solvent=solvent) ", "")
                content = content.replace("max_SCF_cycles", gaussian_parameters["max_SCF_cycles"])
                if dir_parameters[a][0] != "optimization" and dir_parameters[a][0] != "OptRot":
                    content = content.replace("SPECTROSCOPY", dir_parameters[a][0])
                elif dir_parameters[a][0] == "optimization":
                    content = content.replace(" Freq=SPECTROSCOPY", "")
                elif dir_parameters[a][0] == "OptRot":
                    content = content.replace("Freq=SPECTROSCOPY", "Polar=OptRot CPHF=RdFreq")
                if len(dir_parameters[a][3]) != 2 and dir_parameters[a][3] not in custom and dir_parameters[a][3] not in custom1:
                    content = content.replace("BASIS_SET", dir_parameters[a][3])
                else:
                    content = content.replace("BASIS_SET", 'Gen')
                if gaussian_parameters["Overlay"]:
                    content = content.replace("IOp", gaussian_parameters["IOp"])
                else:
                    content = content.replace("IOp", "") 

                file.seek(0)
                file.write(content)
                file.truncate()

                # Writing coordinate data from trajectory, input geometry or Gaussian output.
                for i in range(len(dir_data[a][conformer_count-1][0])):
                    file.write("{:<2} \t {:.6f} \t {:.6f} \t {:.6f}\n".format(dir_data[a][conformer_count-1][1][i], dir_data[a][conformer_count-1][2][i], dir_data[a][conformer_count-1][3][i], dir_data[a][conformer_count-1][4][i]))
                file.write("\n")

                # Writing frequency of incident radiation for ROA.
                if dir_parameters[a][0] == 'ROA' or dir_parameters[a][0] == 'OptRot':
                    for j in range(len(dir_parameters[a][6])):
                        file.write(f"{dir_parameters[a][6][j]} nm ")
                        #if len(dir_parameters[a][6]) != j + 1:
                        #    file.write(", ")
                        #else:
                        #file.write(" nm ")
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
                        params = {'elements': all_atom_types}
                        custom_basis_data = requests.get(base_url + '/api/basis/' + dir_parameters[a][3] + '/format/gaussian94', params=params, headers=headers)
                        basis_data = custom_basis_data.text.split('!----------------------------------------------------------------------\n\n\n')
                        file.write(basis_data[1])
                        if custom_basis_data.status_code != 200:
                            raise RuntimeError("Could not obtain data from the BSE. Check the error information above.")
                        basis = current_basis
                    elif current_basis == basis:
                        file.write(basis_data[1])

                # Pulling custom basis set data from the SoAPy basis set library.
                if dir_parameters[a][3] in custom1 and len(dir_parameters[a][3]) != 2:
                    current_basis = dir_parameters[a][3]
                    if current_basis != basis:
                        basis_file = os.path.join(data_locations, f"basis_sets/{dir_parameters[a][3]}.gbs")
                        basis_file = open(f"{basis_file}", "r")
                        lines = basis_file.readlines()
                        lines = np.asarray(lines)
                        for atom in all_atom_types:
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

                    # Write data for solute basis set.
                    file.write("1 - {} 0\n".format(num_solute_atoms))

                    # Pulling custom basis set data from the Basis Set Exchange API.
                    if dir_parameters[a][3][0] in custom:
                        current_basis = dir_parameters[a][3][0]
                        if current_basis != basis:
                            main_bse_url = "http://basissetexchange.org"
                            base_url = os.environ.get('BSE_API_URL', main_bse_url)
                            headers = {
                                'User-Agent': 'SOlvation Algorithm in PYthon (SoAPy)',
                                'From': 'bshumberger@vt.edu'
                            }
                            params = {'elements': solute_atom_types}
                            custom_basis_data = requests.get(base_url + '/api/basis/' + dir_parameters[a][3][0] + '/format/gaussian94', params=params, headers=headers)
                            basis_data = custom_basis_data.text.split('!----------------------------------------------------------------------\n\n\n')
                            file.write(basis_data[1])
                            if custom_basis_data.status_code != 200:
                                raise RuntimeError("Could not obtain data from the BSE. Check the error information above.")
                            basis = current_basis
                        elif current_basis == basis:
                            file.write(basis_data[1])
                    
                    # Pulling custom basis set data from the SoAPy basis set library.
                    elif dir_parameters[a][3][0] in custom1:
                        current_basis = dir_parameters[a][3][0]
                        if current_basis != basis:
                            basis_file = os.path.join(data_locations, f"basis_sets/{dir_parameters[a][3][0]}.gbs")
                            basis_file = open(f"{basis_file}", "r")
                            lines = basis_file.readlines()
                            lines = np.asarray(lines)
                            for atom in solute_atom_types:
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
                        elif current_basis == basis:
                            file.write(basis_data[1])
                    else:
                        file.write("{}\n".format(dir_parameters[a][3][0]))
                        file.write("****\n")

                    # Write data for solvent basis set.
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
                            params = {'elements': solvent_atom_types}
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

                    # Pulling custom basis set data from the SoAPy basis set library.
                    elif dir_parameters[a][3][1] in custom1:
                        current_basis = dir_parameters[a][3][1]
                        if current_basis != basis:
                            basis_file = os.path.join(data_locations, f"basis_sets/{dir_parameters[a][3][1]}.gbs")
                            basis_file = open(f"{basis_file}", "r")
                            lines = basis_file.readlines()
                            lines = np.asarray(lines)
                            for atom in solvent_atom_types:
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
                        elif current_basis == basis:
                            file.write(basis_data[1])
                            file.write("\n")
                    else:
                        file.write("{}\n".format(dir_parameters[a][3][1]))
                        file.write("****\n")
                        file.write("\n")

                file.write("\n")
                file.write("\n")
                file.write("\n")






