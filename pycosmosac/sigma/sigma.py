import simplejson as json
import numpy as np
from scipy.spatial import distance
from pycosmosac.utils import constants

'''
References:
1. https://dx.doi.org/10.1021/acs.jctc.9b01016
'''


SIGMA_TEMPLATE = """# meta: {meta:s}
# Rows are given as: sigma [e/A^2] followed by a space, then psigmaA [A^2]
# In the case of three sigma profiles, the order is NHB, OH, then OT
"""


def average_charge_density(segments, parameters):
    #Eq. 4 in Ref. 1.
    xyz = segments["xyz"] * constants.au2angstrom 
    charge = segments["charge"]
    area = segments["area"]
    sigma = segments["sigma"]
    if not isinstance(sigma, np.ndarray):
        sigma = charge / area

    rn2 = area / np.pi
    if "r_av" in parameters:
        reff2 = parameters["r_av"]**2
    else:
        reff2 = parameters["a_eff"] / np.pi
    f_decay = parameters["f_decay"]

    d2 = distance.cdist(xyz, xyz)**2
    tmp = np.exp(-f_decay * d2 / (rn2 + reff2)) * rn2 * reff2 / (rn2 + reff2)
    sigma_average = np.sum(tmp * sigma, axis=1) / np.sum(tmp, axis=1)
    return sigma_average

def compute_sigma_profile(segments, sigma, bounds=[-0.025, 0.025], nbins = 50, grid=None):
    #Eqs. 5--8 in Ref. 1.

    step_size = (bounds[1] - bounds[0]) / nbins
    if grid is None:
        grid = np.arange(bounds[0], bounds[1]+step_size, step_size)
    pA = np.zeros((nbins+1), dtype=float)

    area = segments["area"]
    index = (np.floor((sigma - bounds[0]) / step_size)).astype(int)
    w = (grid[index+1] - sigma) / step_size
    for i in range(len(area)):
        pA[index[i]] += w[i] * area[i]
        pA[index[i]+1] += (1 - w[i]) * area[i]
    assert(abs(np.sum(pA) - np.sum(area)) < step_size)
    return pA

class Sigma():
    '''
    Class for sigma profiles
    Attributes:
        mol : Mole
            Molecular information.
        parameters : Parameters
            Parameter information.
        bounds : list
            Boundaries of sigma profile. Default is [-0.025, 0.025].
        nbins : int
            Number of bins in sigma profile. Default is 50.
        write_sigma_file : bool
            Wether or not to write the sigma file. Default is True.
        sigma_filename: str
            Sigma file name. Default is "out.sigma".
        sigma : ndarray
            Averaged charge density.
        pA : ndarray
            y axis of the sigma profile.
        sigma_grid: ndarray
            x axis of the sigma profile.
    '''
    def __init__(self, mol, parameters):
        self.mol = mol
        self.parameters = parameters
        self.bounds = [-0.025, 0.025]
        self.nbins = 50
        self.write_sigma_file = True
        self.sigma_filename = "out.sigma"
        #followings are saved data
        self.sigma = None
        self.pA = None
        self.sigma_grid = None
        self.meta = {}

    def kernel(self):
        self.sigma = average_charge_density(self.mol.cavity.segments, self.parameters.parameters)
        if not (np.amax(self.sigma) <= self.bounds[1] and np.amin(self.sigma) >= self.bounds[0]):
            raise RuntimeError("bounds too narrow. sigma in range [%s, %s]" % (np.amin(self.sigma), np.amax(self.sigma)))

        step_size = (self.bounds[1] - self.bounds[0]) / self.nbins
        self.sigma_grid = np.arange(self.bounds[0], self.bounds[1]+step_size, step_size)
        self.pA = compute_sigma_profile(self.mol.cavity.segments, self.sigma, self.bounds, self.nbins, self.sigma_grid)

        if self.write_sigma_file:
            self.dump_to_file(self.sigma_filename)


    def dump_to_file(self, name):
        '''
        Write the sigma file
        '''
        try:
            with open(name, 'w') as f:
                out = SIGMA_TEMPLATE.format(meta = json.dumps(self.meta, ignore_nan=True))
                for i in range(len(self.pA)):
                    out += '{0:0.4f} {1:17.14e}\n'.format(self.sigma_grid[i], self.pA[i])
                f.write(out)
        except:
            raise RuntimeError("failed to write sigma file.")

if __name__ == "__main__":
    from pycosmosac.param import parameters
    from pycosmosac.cosmo import cosmo

    myparam = parameters.Parameters()
    mycosmo = cosmo.Cosmo()
    mol = mycosmo.load("h2o.cosmo")
    mysigma = Sigma(mol, myparam)
    mysigma.kernel()
