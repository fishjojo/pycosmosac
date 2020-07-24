import numpy as np
from pycosmosac.utils import constants as const 

def calc_G_binary(lngamma, T, vaporP=None, P=None, density=None, molar_mass=None):
    '''
    Compute free energy for binary mixtures from infinite dilution activity coefficients
    '''
    RT = const.R_SI * T
    frac = RT / 1000.0 * const.kj2kcal
    G = frac * lngamma
    if vaporP and density and molar_mass:
        #index 0 for solute; index 1 for solvent
        molarity = np.asarray(density) * 1000.0 / np.asarray(molar_mass)
        G += frac * np.log(vaporP / (RT * molarity[-1]))
        if P and len(density) == 2:
            exp = np.exp((vaporP - P) / (RT * molarity[0]))
            G -= frac * np.log(exp)
    return G
