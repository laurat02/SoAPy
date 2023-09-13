"""Provide functions for analyzing the resultant data."""

import numpy as np
import time 
import os
import shutil
import math
import matplotlib as mpl 
from scipy.signal import find_peaks
import matplotlib.pyplot as plt 



def single_denominator_overlap(dir_frequency_axis, dir_intensity_axis, sample_index, reference_index):
    """
    Function for computing the overlap between two spectra. This function includes a normalization associated with only one of the two spectra rather than both.
    This is specifically designed to include discrepancies between intensities which the "double_denominator_overlap" would neglect.
    """
    # Set sample variables.
    sample_intensity = np.array(dir_intensity_axis[sample_index])

    # Set reference variables.
    reference_frequency = np.array(dir_frequency_axis[reference_index])
    reference_intensity = np.array(dir_intensity_axis[reference_index])

    # Calculate numerator.
    numerator = np.trapz(reference_intensity * sample_intensity, x = reference_frequency)

    # Calculate denominator.
    #if np.trapz(reference_intensity * reference_intensity, x = reference_frequency) > np.trapz(sample_intensity * sample_intensity, x = reference_frequency):
    denominator = np.trapz(reference_intensity * reference_intensity, x = reference_frequency)
    #else:
    #    denominator = np.trapz(sample_intensity * sample_intensity, x = reference_frequency)

    # Calculate overlap.
    overlap = numerator / denominator

    return overlap



def double_denominator_overlap(dir_frequency_axis, dir_intensity_axis, sample_index, reference_index):
    """ 
    Function for computing the overlap between two spectra. This function includes a normalization associated with only one of the two spectra rather than both.
    This is specifically designed to include discrepancies between intensities which the "double_denominator_overlap" would neglect.
    """
    # Set sample variables.
    sample_intensity = np.array(dir_intensity_axis[sample_index])

    # Set reference variables.
    reference_frequency = np.array(dir_frequency_axis[reference_index])
    reference_intensity = np.array(dir_intensity_axis[reference_index])

    # Calculate numerator.
    numerator = np.trapz(reference_intensity * sample_intensity, x = reference_frequency)

    # Calculate denominator.
    denominator = np.sqrt(np.trapz(reference_intensity * reference_intensity, x = reference_frequency) * np.trapz(sample_intensity * sample_intensity, x = reference_frequency))

    # Calculate overlap.
    overlap = numerator / denominator

    return overlap



def integrated_difference(dir_frequency_axis, dir_intensity_axis, sample_index, reference_index):
    """
    Integrates the difference between two functions and makes a ratio of them of the combined area of each spectrum.
    """
    # Set sample variables.
    sample_intensity = np.array(dir_intensity_axis[sample_index])

    # Set reference variables.
    reference_frequency = np.array(dir_frequency_axis[reference_index])
    reference_intensity = np.array(dir_intensity_axis[reference_index])

    # Calculate numerator.
    numerator = np.trapz(reference_intensity**2, x = reference_frequency) - np.trapz(sample_intensity**2, x = reference_frequency)

    # Calculate denominator.
    denominator = np.trapz(reference_intensity**2, x = reference_frequency)

    # Calculate integrated difference.
    int_diff = numerator / denominator

    return int_diff



def compute_statistics(dir_intensities, sample_index, reference_index):
    """
    Computes the differences in mean, variance, and standard deviation with respect to some reference.
    """
    # Calculate mean intensity.
    sample_mean = np.average(dir_intensities[sample_index])
    reference_mean = np.average(dir_intensities[reference_index])

    # Calculate variance in intensities.
    sample_variance_numerator = 0
    reference_variance_numerator = 0
    for b in range(len(dir_intensities[sample_index])):
        sample_variance_numerator += (dir_intensities[sample_index][b] - sample_mean)**2
        sample_variance = sample_variance_numerator / len(dir_intensities[sample_index])

        reference_variance_numerator += (dir_intensities[reference_index][b] - reference_mean)**2
        reference_variance = reference_variance_numerator / len(dir_intensities[reference_index])

    # Calculate the standard deviation.
    sample_standard_deviation = np.sqrt(sample_variance)
    reference_standard_deviation = np.sqrt(reference_variance)

    # Calculate absolute differences.
    mean_diff = abs(reference_mean - sample_mean)
    variance_diff = abs(reference_variance - sample_variance)
    standard_deviation_diff = abs(reference_standard_deviation - sample_standard_deviation)

    return mean_diff, variance_diff, standard_deviation_diff



def wrong_signs(dir_frequencies, dir_intensities, sample_index, reference_index):
    """
    Prints out the frequencies at which there is a sign change. This could be due to normal mode rearrangements or wrong signs. Further developement is needed to distinguish.
    """
    # Set sample variables.
    sample_frequency = np.array(dir_frequencies[sample_index])
    sample_intensity = np.array(dir_intensities[sample_index])

    # Set reference variables.
    reference_frequency = np.array(dir_frequencies[reference_index])
    reference_intensity = np.array(dir_intensities[reference_index])

    # Calculate frequency displacement and sign change.
    delta_freq = []
    sign_change = []
    for a in range(len(reference_frequency)):
        frequency_displacement = reference_frequency[a] - sample_frequency[a]
        intensity_sign = sample_intensity[a] / reference_intensity[a]
        if intensity_sign > 0:
            sign = 1
        elif intensity_sign < 0:
            sign = -1

        delta_freq.append(frequency_displacement)
        sign_change.append(sign)

    return delta_freq, sign_change









