import numpy as np
from pycosmosac.utils import constants as const 

def calc_G_binary(lngamma, T, vaporP=None, P=None, density=None, molar_mass=None):
    '''
    Compute free energy for binary mixtures from infinite dilution activity coefficients
    '''
    RT = const.R_SI * T
    frac = RT / 1000.0 * const.kj2kcal
    G = frac * lngamma
    if vaporP is not None and vaporP < 1e-6 and density and molar_mass:
        molarity = np.asarray(density) * 1000.0 / np.asarray(molar_mass)
        G += frac * np.log(1.0 / (RT * molarity[-1]))
    elif vaporP and density and molar_mass:
        #index 0 for solute; index 1 for solvent
        molarity = np.asarray(density) * 1000.0 / np.asarray(molar_mass)
        G += frac * np.log(vaporP / (RT * molarity[-1]))
        if P and len(density) == 2:
            exp = np.exp((vaporP - P) / (RT * molarity[0]))
            G -= frac * np.log(exp)
    return G

def calc_G_RS(mu_s, mu_ig, T, P, density, molar_mass):
    G = mu_s - mu_ig
    RT = const.R_SI * T
    frac = RT / 1000.0 * const.kj2kcal
    molarity = density * 1000.0 / molar_mass #mol/m^3
    #G -= frac * np.log(molarity * RT / P)
    G -= frac * np.log(molarity * RT)
    return G

if __name__ == "__main__":
    lngamma = 3.0
    T = 300
    G = calc_G_binary(lngamma, T)
    print(G - 1.7884858072629)

    density = [998.]
    molar_mass = [18.]
    vaporP = 1e4
    G = calc_G_binary(lngamma, T, vaporP=vaporP, density=density, molar_mass=molar_mass)
    print(G - -3.8956651094414827)

    density = [800., 998.]
    molar_mass = [20., 18.]
    vaporP = 1e4
    P = 1e5
    G = calc_G_binary(lngamma, T, vaporP=vaporP, P=P, density=density, molar_mass=molar_mass)
    print(G - -3.8951273459414826)

