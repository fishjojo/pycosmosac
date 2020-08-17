import os
import numpy as np
from pycosmosac.param import data, parameters
from pycosmosac.cosmo import cosmo
from pycosmosac.sigma import sigma
from pycosmosac.ac import ac
from pycosmosac.utils import thermo

'''
An example to compute hydration free energy of carbofuran
'''

path = os.path.abspath(os.path.dirname(__file__)) + "/"

#load parameters
myparam = parameters.Parameters(parameters = data.BIOSAC_SVP_GEPOL)

#compute sigma profile for solute
solute = "carbofuran"
mol1 = cosmo.Cosmo().load(path + solute + ".cosmo")
sigma1 = sigma.Sigma(mol1, myparam)
#whether or not to write the sigma file
sigma1.write_sigma_file = False
sigma1.sigma_filename = solute + ".sigma"
#default bounds for the grid; if "bounds too narrow" error occures, increase the numbers
sigma1.bounds = [-0.025, 0.025]
sigma1.kernel()

#compute sigma profile for solvent
solvent = "h2o"
mol2 = cosmo.Cosmo().load(path + solvent + ".cosmo")
sigma2 = sigma.Sigma(mol2, myparam)
sigma2.write_sigma_file = False
sigma2.sigma_filename = solvent + ".sigma"
#bounds should be consistent for every component in the solution
sigma2.bounds = [-0.025, 0.025]
sigma2.kernel()

#mole fraction of solute and solvent; [0, 1] means infinite dilution
x = [0.0, 1.0]
#temperature in Kelvin
T = 298.15
#compute activity coefficient of solute
myac = ac.AC([mol1,mol2], x, T, [sigma1,sigma2], myparam)
lngamma = myac.kernel()[0]

#vapor pressure of solute in Pascal
vp_solute = 0.0006466
#density of solute and solvent in g/L
density = [np.nan, 997.0]
#molar mass of solute and solvent in g/mol
molar_mass=[np.nan, 18.0153]

print("experimental hydration free energy: -9.61 kcal/mol")
#compute delta_G in kcal/mol
G = thermo.calc_G_binary(lngamma, T, vaporP = vp_solute, density = density, molar_mass = molar_mass)
print("computed hydration free energy (using experimental vapor pressure): {:.2f} kcal/mol".format(G))

#if vapor pressure is unknown and only the relative free energies of similar compounds are of interest,
#then it's reasonable to set vapor pressure as zero.
G0 = thermo.calc_G_binary(lngamma, T, vaporP = 0.0, density = density, molar_mass = molar_mass)
print("computed hydration free energy (using 0 vapor pressure): {:.2f} kcal/mol".format(G0))

#if vapor pressure is unknown and the absolute hydration energy is of interest, 
#then we can use an empirical model to estimate the contribution from the gas-phase chemical potential.
#note this may introduce a few kcal/mol error, and is only recommended for molecules with small vapor pressures, 
#because those were the only systems considered in the parameter fitting.
E_vp = ac.e_vapor(mol1, sigma1, 0., 0., T, myparam.parameters, data.BIOSAC_SVP_GEPOL_disp, True)
G = G0 + E_vp
print("computed hydration free energy (using estimated vapor pressure): {:.2f} kcal/mol".format(G))
