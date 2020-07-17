import warnings
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

def compute_sigma_profile(area, sigma, bounds=[-0.025, 0.025], nbins = 50, grid=None):
    #Eqs. 5--8 in Ref. 1.

    step_size = (bounds[1] - bounds[0]) / nbins
    if grid is None:
        grid = np.arange(bounds[0], bounds[1]+step_size, step_size)
    pA = np.zeros((nbins+1), dtype=float)

    index = (np.floor((sigma - bounds[0]) / step_size)).astype(int)
    w = (grid[index+1] - sigma) / step_size
    for i in range(len(area)):
        pA[index[i]] += w[i] * area[i]
        pA[index[i]+1] += (1 - w[i]) * area[i]
    assert(abs(np.sum(pA) - np.sum(area)) < step_size)
    return pA

def split_sigma_profile(mol, sigma, parameters, bounds=[-0.025, 0.025], nbins = 50, grid=None):
    step_size = (bounds[1] - bounds[0]) / nbins
    if grid is None:
        grid = np.arange(bounds[0], bounds[1]+step_size, step_size)

    geometry = mol.geometry
    atoms = np.array(geometry["atom"])
    hb_class = np.array(mol.hb_class)
    atom_map = mol.cavity.atom_map
    segments = mol.cavity.segments
    area = segments["area"]

    seg_atoms = atoms[atom_map-1]
    seg_hb_class = hb_class[atom_map-1]

    mask_oh = (((seg_atoms == "O") & (seg_hb_class == "OH") & (sigma > 0.0)) |
               ((seg_atoms == "H") & (seg_hb_class == "OH") & (sigma < 0.0)))
    mask_ot = (((np.isin(seg_atoms, ["O","N","F"])) & (seg_hb_class == "OT") & (sigma > 0.0)) |
               ((seg_atoms == "H") & (seg_hb_class == "OT") & (sigma < 0.0)))
    mask_nhb = ~(mask_oh | mask_ot)

    sigma_oh = sigma[mask_oh]
    sigma_ot = sigma[mask_ot]
    sigma_nhb = sigma[mask_nhb]
    area_oh = area[mask_oh]
    area_ot = area[mask_ot]
    area_nhb = area[mask_nhb]

    pA_oh = compute_sigma_profile(area_oh, sigma_oh, bounds, nbins, grid)
    pA_ot = compute_sigma_profile(area_ot, sigma_ot, bounds, nbins, grid)
    pA_nhb = compute_sigma_profile(area_nhb, sigma_nhb, bounds, nbins, grid)

    if not "sigma_0" in parameters:
        warings.warn("sigma_0 missing in the parameter set. probabality of hydrogen bonding will be set to 1")
        P_hb=1.0
    else:
        sigma0 = parameters["sigma_0"]
        P_hb = 1.0 - np.exp(-grid*grid/(2.0*sigma0*sigma0))

    pA_hb = pA_oh + pA_ot
    pA_oh *= P_hb
    pA_ot *= P_hb
    pA_nhb = pA_nhb + pA_hb * (1.0 - P_hb)
    return pA_nhb, pA_oh, pA_ot


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
        sigma_filename : str
            Sigma file name. Default is "out.sigma".
        split_sigma : bool
            Wether to split sigma profile. Default is False.
        chem_name : str
            Name of the compound.
        cas_no : str
            CAS registry number.
        inchikey : str
            IUPAC Standard InChIKey
        sigma : ndarray
            Averaged charge density.
        pA : ndarray
            y axis of the sigma profile.
        pA_nhb : ndarray
            Non-hydrogen bonding part of sigma profile.
        pA_oh : ndarray
            OH type of hydrogen bonding part of sigma profile.
        pA_ot : ndarray
            Other type of hydrogen bonding part of sigma profile.
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
        self.split_sigma = False
        #optional inputs
        self.chem_name = None
        self.cas_no = None
        self.inchikey = None
        #followings are saved data
        self.sigma = None
        self.pA = None
        self.pA_nhb = None
        self.pA_oh = None
        self.pA_ot = None
        self.sigma_grid = None
        self.meta = {}

    def kernel(self):
        self.sigma = average_charge_density(self.mol.cavity.segments, self.parameters.parameters)
        if not (np.amax(self.sigma) <= self.bounds[1] and np.amin(self.sigma) >= self.bounds[0]):
            raise RuntimeError("bounds too narrow. sigma in range [%s, %s]" % (np.amin(self.sigma), np.amax(self.sigma)))

        step_size = (self.bounds[1] - self.bounds[0]) / self.nbins
        self.sigma_grid = np.arange(self.bounds[0], self.bounds[1]+step_size, step_size)
        self.pA = compute_sigma_profile(self.mol.cavity.segments["area"], self.sigma, self.bounds, self.nbins, self.sigma_grid)
        if self.split_sigma:
            self.pA_nhb, self.pA_oh, self.pA_ot = split_sigma_profile(self.mol, self.sigma, self.parameters.parameters,
                                                                      self.bounds, self.nbins, self.sigma_grid)

        if self.write_sigma_file:
            self.dump_to_file(self.sigma_filename)

    def get_meta(self):
        meta = {}
        if self.chem_name:
            meta['name'] = self.chem_name
        if self.cas_no:
            meta['CAS'] = self.cas_no
        if self.inchikey:
            meta["standard_INCHIKEY"] = self.inchikey
        meta['area [A^2]'] = self.mol.cavity.area
        meta['volume [A^3]'] = self.mol.cavity.volume

        parameters = self.parameters.parameters
        if "r_av" in parameters:
            meta['r_av [A]'] = parameters["r_av"]
        else:
            meta['r_av [A]'] = (parameters["a_eff"] / np.pi)**0.5
        meta['f_decay'] = parameters["f_decay"]
        if "sigma_hb" in parameters:
            meta['sigma_hb [e/A^2]'] = parameters["sigma_hb"]
        return meta

    def dump_to_file(self, name):
        '''
        Write the sigma file
        '''
        try:
            with open(name, 'w') as f:
                out = SIGMA_TEMPLATE.format(meta = json.dumps(self.get_meta(), ignore_nan=True))
                if self.split_sigma:
                    for profiles in [self.pA_nhb, self.pA_oh, self.pA_ot]:
                        for i in range(len(profiles)):
                            out += '{0:0.4f} {1:17.14e}\n'.format(self.sigma_grid[i], profiles[i])
                else:
                    for i in range(len(self.pA)):
                        out += '{0:0.4f} {1:17.14e}\n'.format(self.sigma_grid[i], self.pA[i])
                f.write(out)
        except:
            raise RuntimeError("failed to write sigma file.")

if __name__ == "__main__":
    from pycosmosac.param import parameters, data
    from pycosmosac.cosmo import cosmo
    from pycosmosac.utils import misc

    mycosmo = cosmo.Cosmo()
    mol = mycosmo.load("./test/h2o.cosmo")
    myparam = parameters.Parameters()
    mysigma = Sigma(mol, myparam)
    mysigma.write_sigma_file = False
    mysigma.kernel()
    print(misc.fp(mysigma.pA) - -1.917622937152116)

    myparam = parameters.Parameters(data.Hsieh_2010)
    mysigma = Sigma(mol, myparam)
    mysigma.write_sigma_file = False
    mysigma.split_sigma = True
    mysigma.kernel()
    print(np.sum(mysigma.pA - (mysigma.pA_nhb + mysigma.pA_oh + mysigma.pA_ot)))
