# A from scratch Hartree-Fock program
# Imports
import numpy as np
from scipy.special import erf
from hf_file import xyz_reader

xyz_reader('HeH.xyz')
