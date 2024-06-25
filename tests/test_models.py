import numpy as np
import pytest
import sys
import os

#TODO: fix this import 
loc = '/Users/allybaldelli/Desktop/AMNH-stuff/Doublet-work/Doublet-Quantifiers'
sys.path.append(os.path.abspath(loc))
from Doublet_Quantifiers.models import *

def test_gaussian():
    """Test the gaussian function outputs the correct values."""
    x = np.arange(0, 10, 2)
    A = 1
    sigma = 1
    mu = 0
    output_array = gaussian(x, A, sigma, mu)
    expected_array = np.array([1.00000000e+00, 1.35335283e-01, 3.35462628e-04, 1.52299797e-08,
       1.26641655e-14])
    np.testing.assert_allclose(output_array, expected_array)

# def test_voigt():
#     """Test the voigt function outputs the correct values."""
#     x = np.arange(0, 10, 2)
#     A = 1
#     mu = 0
#     sigma = 1
#     gamma = 1
#     output_array = voigt(x, A, sigma, gamma, mu)
#     expected_array = np.array([])
#     np.testing.assert_allclose(output_array, expected_array)