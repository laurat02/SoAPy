import matplotlib.pyplot as plt
from SoAPy.functions import set_options
from SoAPy.functions import convert_GROMACS
from SoAPy.functions import generate_files
from SoAPy.functions import generate_coordinates
from SoAPy.functions import generate_batch_submission_script
from SoAPy.functions import collect_data
from SoAPy.functions import generate_spectrum
import time
import os

cwd = os.getcwd()

# Location of MD trajectory.
trajectory_location = cwd + '/' + 'output_data.txt'

# Molecule name
molecule_name = 'B_glucose'

# Set testing parameters in lists.
spectroscopy = ['VCD']
shell_type = ['molecular']
functional = ['CAM-B3LYP']
basis_set = ['6-31G(d)', 'cc-pVDZ', ['cc-pVDZ', '6-31G(d)']]
distance_threshold = ['2.5']
snapshots = ['10']
frequency = []

# Set spectrum parameters.
fwhm = 0.001
number_of_points = 2000

t0 = time.time()

# What type of test is being researched, basis sets, distance thresholds, or number of snapshots?
dir_list, dir_parameters, relative_dir_list = set_options(spectroscopy, shell_type, functional, basis_set, distance_threshold, snapshots, frequency)
t1 = time.time()

# Collect frequency and intensity data from Gaussian output files.
dir_frequencies, dir_intensities = collect_data(cwd, dir_list, dir_parameters)
t2 = time.time()

# Generate spectral data.
dir_frequency_axis, dir_intensity_axis = generate_spectrum(fwhm, number_of_points, dir_list, dir_parameters, dir_frequencies, dir_intensities)
t3 = time.time()

# Plot spectra.
# Setting plot parameters.
plt.rcParams["figure.figsize"] = [13.50, 8.50]
plt.rcParams["figure.autolayout"] = True

# Selecting data to plot.
for a in range(len(dir_list)):
    plt.plot(dir_frequency_axis[a], dir_intensity_axis[a], label=f"{dir_parameters[a][3]}")

    # Plotting the spectrum with x- and y- axis labels.
    plt.xlabel('Frequency (cm-1)')

    if dir_parameters[a][0] == 'VCD':
        plt.ylabel('Differential Molar Extinction Coefficient')
    if dir_parameters[a][0] == 'ROA':
        plt.ylabel('ROA Intensity Difference')

#Plotting legend.
plt.legend(loc="upper left")

# Setting the limits on the x-axis.
plt.xlim(500,1100)

# Saving the plot.
plt.savefig(f'/Users/bshumberger/Documents/SOlvation_Algorithm/spectrum.pdf', bbox_inches='tight')

print("Total Time: ", t3-t0)








