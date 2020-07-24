import numpy as np
from pycosmosac.utils import constants as const 

def calc_G_binary(lngamma, vaporP, P, T, solvent_density, solvent_molar_mass):
    molarity = solvent_density * 1000.0 / solvent_molar_mass
    RT = const.R_SI * T
    frac = RT / 1000.0 * const.kj2kcal
    exp = np.exp((vaporP - P) / (RT * molarity))
    G = (lngamma + np.log(vaporP / (RT * molarity * exp))) * frac
    return G

