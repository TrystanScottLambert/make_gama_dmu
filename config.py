"""
Configuration settings in python module
"""

from nessie import FlatCosmology
from astropy.cosmology import FlatLambdaCDM
import numpy as np

cosmo = FlatCosmology(0.7, 0.3)
astropy_cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
OVERFACTOR = 400  # The factor the randoms were expanded.
B0, R0 = 0.06, 32
APPARENT_MAG_LIM = 19.65
AB_CUT = -10
MASS_A = 10
MASS_FUNC_PARAMS = np.array([2.0, 17.9, 1.5])
LUM_B = 1.04
LUM_FUNC_PARAMS = np.array([0.65, -0.5, 0.22])
SUN_MAG = 4.67

GROUP_ID_OFFSET = {"g09": 1e5, "g12": 2e5, "g15": 3e5, "g23": 5e5}

NAME = "Trystan Lambert"
EMAIL = "trystan.lambert@uwa.edu.au"
VERSION = "v11"
ABS = "Galaxy groups for the G09, G12, G15, and G23 regions.\n# Groups built using Nessie and new GAMA III photometry based on profound (gkInputCat)"