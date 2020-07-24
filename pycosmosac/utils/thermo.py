import numpy as np
from pycosmosac.utils import constants as const 

def calc_G_binary(lngamma, vaporP, P, T, density, molar_mass):
    molarity = np.asarray(density) * 1000.0 / np.asarray(molar_mass)
    RT = const.R_SI * T
    frac = RT / 1000.0 * const.kj2kcal
    exp = np.exp((vaporP - P) / (RT * molarity[0]))
    G = (lngamma + np.log(vaporP / (RT * molarity[1] * exp))) * frac
    return G

