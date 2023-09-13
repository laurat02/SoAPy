"""Provide functions for plotting the analized data."""

import numpy as np
import time 
import os
import shutil
import math
import matplotlib as mpl 
import matplotlib.pyplot as plt 



def num_tests(dir_parameters, compare='basis_set'):
    if compare == 'basis_set':
        spec_list = []
        func_list = []

        # Get number of spectroscopies and density functionals.
        for a in range(len(dir_parameters)):
            spec_list.append(dir_parameters[a][0])
            func_list.append(dir_parameters[a][2])
        spec_unique = set(spec_list)
        func_unique = set(func_list)
        spec = list(spec_unique)
        func = list(func_unique)
    
        num_tests = len(spec) * len(func)

    else:
        print("More tests comparisons coming soon.")

    return num_tests, spec, func

def normalize_by_test(dir_parameters, compare='basis_set'):
    number_of_tests, spectroscopies, functionals = num_tests(dir_parameters)
    
    # Setting the beggining and ending indices for each test.
    for a in range(1, number_of_tests + 1):
        begin = a * len(basis_set) - len(basis_set)
        end = a * len(basis_set)

def plot_spectra(compare = 'basis_set', )
