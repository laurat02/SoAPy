"""Provide the primary functions."""

import numpy as np
import time
import os
import shutil
import math
import requests
import matplotlib as mpl
import matplotlib.pyplot as plt
from statistics import mean


def collect_VCD_and_ROA_data(cwd, dir_list, dir_parameters):
    """
    Collects the VCD and ROA data from the Gaussian output files for spectral generation.
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
            test_frequencies.append(real_frequencies)
            test_intensities.append(real_intensities)

            conformer_count += 1

        dir_frequencies.append(test_frequencies)
        dir_intensities.append(test_intensities)
        dir_nbf.append(nbf)

    return dir_frequencies, dir_intensities, dir_nbf



def collect_OptRot_data(cwd, dir_list, dir_parameters):
    """
    Collects the data optical rotation data from the Gaussian output files.
    """
    # Initialize directory level arrays.
    dir_frequencies = []
    dir_intensities = []
    dir_molar_mass = []
    dir_nbf = []

    # Change to test specific directory.
    for a in range(len(dir_list)):
        print("____________________________")
        print(f"Spectroscopy: {dir_parameters[a][0]} || Solvent Shell Type: {dir_parameters[a][1]} || Functional: {dir_parameters[a][2]} || Basis: {dir_parameters[a][3]} || Distance: {dir_parameters[a][4]} || Snapshots: {dir_parameters[a][5]}")
        with open(f"{cwd}/output_data.txt", "a") as file:
            file.write(f"Spectroscopy: {dir_parameters[a][0]} || Solvent Shell Type: {dir_parameters[a][1]} || Functional: {dir_parameters[a][2]} || Basis: {dir_parameters[a][3]} || Distance: {dir_parameters[a][4]} || Snapshots: {dir_parameters[a][5]}\n")

        test_frequencies = []
        test_intensities = []
        test_molar_mass = []

        # Obtain data from each conformer in the test.
        conformer_count = 1
        while conformer_count <= int(dir_parameters[a][5]):
            #print(f"Snapshot: {conformer_count}")
            frequency = []
            intensity = []
            molar_mass = []
            os.chdir(f"{dir_list[a]}/cmpd_{conformer_count}")

            # Data collection for optical rotation.
            with open("output.log", "r") as file_out:
                for line in file_out:
                    split_line = line.split()
                    if len(split_line) > 7:
                        if split_line[0] == "NAtoms=":
                            natom = int(split_line[1])
                        if split_line[0] == "Molar":
                            frequency.append(float(split_line[7]))
                            intensity.append(float(split_line[10]))
                            molar_mass.append(float(split_line[3]))
                        if split_line[0] == "NBasis=":
                            nbf = int(split_line[1])

                # Correcting units on incident radiation back to nanometers.
                frequency = [i / 10 for i in frequency]

            print(f"Number of Atoms = {natom} \t Number of Basis Functions = {nbf}")

            # Print frequencies and intensities to output file.
            with open(f"{cwd}/output_data.txt", "a") as file:
                file.write(f"Snapshot = {conformer_count} \t Number of Atoms = {natom}\n")
                for i in range(len(frequency)):
                    file.write("{:.4f} \t {:.4f} \t {:e}\n".format(molar_mass[i], frequency[i], intensity[i]))

            # Appending to test level arrays.
            test_frequencies.append(frequency)
            test_intensities.append(intensity)
            test_molar_mass.append(molar_mass)

            conformer_count += 1

        dir_frequencies.append(test_frequencies)
        dir_intensities.append(test_intensities)
        dir_molar_mass.append(test_molar_mass)
        dir_nbf.append(nbf)

    return dir_frequencies, dir_intensities, dir_molar_mass, dir_nbf



def generate_VCD_and_ROA_spectra(fwhm, number_of_points, dir_list, dir_parameters, dir_frequencies, dir_intensities, min_frequency, max_frequency):
    """
    Generates the spectrum for each test. The spectral generation uses a simple averaging algorithm to weight the snapshots.
    """
    # Initialize directory level arrays for plotting.
    dir_frequency_axis = []
    dir_intensity_axis = []
    max_intensity = []
    dir_frequencies_new = []
    dir_intensities_new = []

    for a in range(len(dir_list)):
        dir_frequencies_new.append([])
        dir_intensities_new.append([])
        for b in range(0, len(dir_frequencies[a])):
            dir_frequencies_new[a].extend(dir_frequencies[a][b])
            dir_intensities_new[a].extend(dir_intensities[a][b])

        # Sorts the frequencies and intensities for a given test in ascending order of frequencies.
        spec_frequencies, spec_intensities = zip(*sorted(zip(dir_frequencies_new[a], dir_intensities_new[a])))

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



def generate_VCD_and_ROA_convergence_spectra(fwhm, number_of_points, dir_list, dir_parameters, dir_frequencies, dir_intensities, min_frequency, max_frequency, snapshot_list):
    """
    Generates the spectrum for each test assuming only one snapshot test has been performed and only VCD and ROA are being explored (not optical rotation). 
    This function returns the running average of data for the chosen interval of snapshots through the total snapshot number. 
    The spectral generation uses a simple averaging algorithm to weight the snapshots.
    """
    # Initialize directory level arrays for plotting.
    dir_frequency_axis = []
    dir_intensity_axis = []
    max_intensity = []
    dir_frequencies_new = []
    dir_intensities_new = []

    for a in range(len(dir_list)):
        print("Generating spectral data for test", a,".")
        dir_frequencies_new.append([])
        dir_intensities_new.append([])
        dir_frequency_axis.append([])
        dir_intensity_axis.append([])
        max_intensity.append([])

        for x in range(0,len(snapshot_list)):
            print("Calculating spectrum for",snapshot_list[x],"snapshots.")
            dir_frequencies_new[a].append([])
            dir_intensities_new[a].append([])

            for b in range(0, snapshot_list[x]):
                dir_frequencies_new[a][x].extend(dir_frequencies[a][b])
                dir_intensities_new[a][x].extend(dir_intensities[a][b])

            #print(dir_frequencies_new[a][x], "\n")
            # Sorts the frequencies and intensities for a given test in ascending order of frequencies.
            spec_frequencies, spec_intensities = zip(*sorted(zip(dir_frequencies_new[a][x], dir_intensities_new[a][x])))

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
                normalized_spec_intensities[j] = spec_intensities[j]*(1/int(snapshot_list[x]))

            # Fitting data to line shapes.
            # See equations 1.16 and 3.51 "Vibrational Optical Activity Principles and Applications" by Laurence Nafie for details.
            # Note the prefactor in equation 3.51 is excluded.
            if dir_parameters[a][0] == 'VCD':
                for c in range(len(frequency_axis)):
                    for d in range(len(spec_frequencies)):
                        # Equation 8d in "ECD Cotton Effects Approximated by the Gaussian Curve and Other Methods" by Philip J. Stephens and Nobuyuki Harada.
                        #intensity_axis[b] += (1/((2.296*10**(-39))*np.sqrt(np.pi)*fwhm))*spec_frequencies_eV[c]*normalized_spec_intensities[c]*np.exp(-((frequency_axis_eV[b]-spec_frequencies_eV[c])/fwhm)**2)
                        intensity_axis[c] += spec_frequencies_eV[d]*normalized_spec_intensities[d]*(fwhm**2/(4*(frequency_axis_eV[c]-spec_frequencies_eV[d])**2+fwhm**2))
            if dir_parameters[a][0] == 'ROA':
                for c in range(len(frequency_axis)):
                    for d in range(len(spec_frequencies)):
                        #intensity_axis[b] += normalized_spec_intensities[c]*(2/np.pi)*(fwhm/(4*(frequency_axis_eV[b]-spec_frequencies_eV[c])**2+fwhm**2))
                        intensity_axis[c] += spec_frequencies_eV[d]*normalized_spec_intensities[d]*(fwhm**2/(4*(frequency_axis_eV[c]-spec_frequencies_eV[d])**2+fwhm**2))

            # NEED TO CHECK UNITS FOR BOTH OF THESE EQUATIONS. THIS INCLUDES UNITS ON FWHM.
            # NEED TO ADD NORMALIZATION FUNCTION FOR THESE EQUATIONS.

            # Appending maximum intensity for normalization.
            max_intensity[a].append(max(np.absolute(intensity_axis)))

            # Appending frequencies and intensities to directory level arrays.
            dir_frequency_axis[a].append(frequency_axis)
            dir_intensity_axis[a].append(intensity_axis)

    return dir_frequency_axis, dir_intensity_axis, max_intensity



def generate_OptRot_convergence_data(dir_list, dir_parameters, dir_frequencies, dir_intensities, dir_molar_mass, snapshot_list, molecule_mass):
    """
    Generates the optical rotation data as a running average running average of data for the chosen interval of snapshots through the total snapshot number.
    Given that various snapshot can have a variable number of molecules, we choose to multiply the optical rotation value, alpha, by the molecular mass.
    """
    # Initialize directory level arrays for plotting.
    dir_frequency_axis = []
    dir_intensity_axis = []
    dir_molar_mass_axis = []
    dir_frequencies_new = []
    dir_intensities_new = []
    dir_snapshot_axis = []
    dir_optrot_axis = []

    # Loop through the tests.
    for a in range(len(dir_list)):
        print("Generating optical rotation data for test", a,".")
        dir_frequencies_new.append([])
        dir_intensities_new.append([])
        dir_frequency_axis.append([])
        dir_intensity_axis.append([])
        dir_snapshot_axis.append([])
        dir_optrot_axis.append([])

        # Loop through the various incident frequencies to create a structure amenable to graphing.
        for e in range(len(dir_parameters[a][6])):
            dir_snapshot_axis[a].append([])
            dir_optrot_axis[a].append([])

        # Loop through the snapshots of interests.
        for x in range(0,len(snapshot_list)):
            print("Calculating averages optical rotation data for",snapshot_list[x],"snapshots.")
            dir_frequencies_new[a].append([])
            dir_intensities_new[a].append([])
            dir_frequency_axis[a].append([])
            dir_intensity_axis[a].append([])

            # Append the previous snapshots to the list for convergence data.
            for b in range(0, snapshot_list[x]):
                dir_frequencies_new[a][x].extend(dir_frequencies[a][b])
                dir_intensities_new[a][x].extend(np.array(dir_intensities[a][b])*dir_molar_mass[a][b][0]/molecule_mass)

            average_frequency = []
            average_intensity = []

            # Loops for averaging across the different frequencies of interest.
            for e in range(len(dir_parameters[a][6])):
                average_frequency.append([])
                average_intensity.append([])
                for d in range(len(dir_frequencies_new[a][x])):
                    if float(dir_frequencies_new[a][x][d]) == float(dir_parameters[a][6][e]):
                        average_frequency[e].append(dir_frequencies_new[a][x][d])
                        average_intensity[e].append(dir_intensities_new[a][x][d])
             
                dir_frequency_axis[a][x].append(mean(average_frequency[e]))
                dir_intensity_axis[a][x].append(mean(average_intensity[e]))

                dir_snapshot_axis[a][e].append(snapshot_list[x])
                dir_optrot_axis[a][e].append(mean(average_intensity[e]))

            #print(dir_frequencies_new[a][x])
            #print(dir_intensities_new[a][x], "\n")

    return dir_frequency_axis, dir_intensity_axis, dir_snapshot_axis, dir_optrot_axis

