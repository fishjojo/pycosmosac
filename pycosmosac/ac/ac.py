import warnings
import numpy as np
from pycosmosac.utils import constants as const

def lngamma_c(mols, x, parameters):
    '''
    Staverman−Guggenheim combinatorial term
    '''
    r0 = parameters["r_0"]
    q0 = parameters["q_0"]
    z = float(parameters["z"])

    A = []
    V = []
    for mol in mols:
        A.append(mol.cavity.area)
        V.append(mol.cavity.volume)
    A = np.asarray(A)
    V = np.asarray(V)

    r = V / r0
    q = A / q0

    xq = x * q
    xr = x * r
    theta_phi = q / r / np.sum(xq) * np.sum(xr)
    phi_x = r / np.sum(xr)
    l = 0.5 * z * (r - q) - (r - 1.0)

    result = np.log(phi_x) + 0.5 * z * q * np.log(theta_phi) + l - phi_x * np.sum(x*l)
    return result

def solve_lnGamma(W, GammaS, ps, T, thresh=1e-3, maxiter=500):
    #TODO improve convergence speed
    count = 0
    while (count < maxiter):
        GammaS_new = np.exp(-np.log(np.sum(np.exp(-W / (const.R*T)) * GammaS * ps, axis=1)))
        GammaS = 0.5 * (GammaS + GammaS_new)
        diff = np.amax(np.absolute(GammaS_new - GammaS))
        if diff < thresh:
            break
        count += 1
    if count == maxiter:
        warnings.warn("solving Gamma reached max iterations: %s steps with diff %.3e" % (maxiter,diff))
    return np.log(GammaS)

def lngamma_r(mols, sigmas, x, T, parameters, thresh=1e-3, maxiter=500):
    aeff = parameters["a_eff"]
    A = []
    for mol in mols:
        A.append(mol.cavity.area)
    A = np.asarray(A)
    n = A / aeff

    p = []
    for i, sigma in enumerate(sigmas):
        p.append(sigma.pA / A[i])

    aa = np.sum(x * A)
    ps = np.zeros_like(p[0])
    for i in range(len(x)):
        ps += x[i] * A[i] * p[i] / aa

    grid = sigmas[0].sigma_grid

    alpha = parameters["alpha_prime"]
    sigma_hb = parameters["sigma_hb"]
    chb = parameters["c_hb"]
    W = 0.5 * alpha * (grid[:, None] + grid[None, :])**2
    sigma_acc = np.maximum(grid[:, None], grid[None, :]) - sigma_hb
    sigma_acc = np.maximum(sigma_acc, 0)
    sigma_don = np.minimum(grid[:, None], grid[None, :]) + sigma_hb
    sigma_don = np.minimum(sigma_don, 0)
    W = W + chb * sigma_acc * sigma_don

    GammaS = np.ones_like(ps)
    lnGammaS = solve_lnGamma(W, GammaS, ps, T, thresh, maxiter)

    result = []
    for pi in p:
        Gammai = np.ones_like(ps)
        lnGammai = solve_lnGamma(W, Gammai, pi, T, thresh, maxiter)
        result.append(np.dot(pi, lnGammaS - lnGammai))

    result = np.asarray(result)
    result = result * n
    return result

def lngamma_r3(mols, sigmas, x, T, parameters, thresh=1e-3, maxiter=500):
    aeff = parameters["a_eff"]
    A = []
    for mol in mols:
        A.append(mol.cavity.area)
    A = np.asarray(A)
    n = A / aeff

    p = []
    for i, sigma in enumerate(sigmas):
        pA = np.asarray([sigma.pA_nhb, sigma.pA_oh, sigma.pA_ot]).ravel()
        p.append(pA / A[i])

    aa = np.sum(x * A)
    ps = np.zeros_like(p[0])
    for i in range(len(x)):
        ps += x[i] * A[i] * p[i] / aa

    try:
        ces = parameters["A_es"] + parameters["B_es"] / (T*T)
        c_ohoh = parameters["c_ohoh"]
        c_otot = parameters["c_otot"]
        c_ohot = parameters["c_ohot"]
    except:
        raise RuntimeError("A_es, B_es, c_ohoh, c_otot, or c_ohot not in parameter set")

    grid = sigmas[0].sigma_grid
    grid3 = np.asarray([grid,grid,grid]).ravel()
    ngrid = len(grid)
    ngrid2 = ngrid * 2
    ngrid3 = len(grid3)
    assert(ngrid3 % ngrid == 0)
    chb = np.zeros((ngrid3,ngrid3))
    sigma2 = grid3[:,None] * grid3[None,:]
    mask_sigma = sigma2 < 0
    mask_oh = np.full((ngrid3,ngrid3), False)
    mask_oh[ngrid:ngrid2, ngrid:ngrid2] = True
    mask_oh = ((mask_oh) & (mask_sigma))
    chb[mask_oh] = c_ohoh

    mask_ot = np.full((ngrid3,ngrid3), False)
    mask_ot[ngrid2:, ngrid2:] = True
    mask_ot = ((mask_ot) & (mask_sigma))
    chb[mask_ot] = c_otot

    mask_ohot = np.full((ngrid3,ngrid3), False)
    mask_ohot[ngrid:ngrid2, ngrid2:] = True
    mask_ohot[ngrid2:, ngrid:ngrid2] = True
    mask_ohot = ((mask_ohot) & (mask_sigma))
    chb[mask_ohot] = c_ohot

    W = ces * (grid3[:, None] + grid3[None, :])**2 - chb * (grid3[:, None] - grid3[None, :])**2

    GammaS = np.ones_like(ps)
    lnGammaS = solve_lnGamma(W, GammaS, ps, T, thresh, maxiter)

    result = []
    for pi in p:
        Gammai = np.ones_like(ps)
        lnGammai = solve_lnGamma(W, Gammai, pi, T, thresh, maxiter)
        result.append(np.dot(pi, lnGammaS - lnGammai))

    result = np.asarray(result)
    result = result * n
    return result

def lngamma_dsp(mols, sigmas, x, parameters):
    if len(x) != 2:
        warnings.warn("dispersion corrections are only implemented for binary mixtures.")
        return 0.0

    sigma_hb = parameters["sigma_hb"]

    grid = sigmas[0].sigma_grid
    disp_epss = []
    disp_types = []
    for i, mol in enumerate(mols):
        disp_eps, disp_type = mol.get_dispersion_type()
        if disp_type == 'NHB':
            mask = (sigmas[i].pA > 0) #I use total profile here
            sigma_max = np.max(grid[mask])
            sigma_min = np.min(grid[mask])
            if sigma_max > sigma_hb:
                if sigma_min > -sigma_hb:
                    disp_type = "HB-ACCEPTOR"
                else:
                    disp_type = "HB-DONOR-ACCEPTOR"
        disp_epss.append(disp_eps)
        disp_types.append(disp_type)

    from pycosmosac.param import data
    w = data.disp["w"]
    A = w * (0.5 * (disp_epss[0]+disp_epss[1]) - (disp_epss[0]*disp_epss[1])**0.5)
    if "H2O" in disp_types and "COOH" in disp_types:
        A *= -1.0
    elif "H2O" in disp_types and "HB-ACCEPTOR" in disp_types:
        A *= -1.0
    elif "COOH" in disp_types and ("HB-DONOR-ACCEPTOR" in disp_types or 'NHB' in disp_types):
        A *= -1.0
    result = np.asarray([A*x[1]**2, A*x[0]**2])
    return result

class AC():
    '''
    Class for computing activity coefficients
    Atrributes:
        mols : list
            List of Mole objects.
        x : ndarray
            Mole fractions.
        T : float
            Temperature in K.
        sigmas : list
            List of Sigma objects.
        parameters : Parameters
            The Parameters object.
        split_sigma : bool
            Whether to use splitted sigma profile.
        dispersion : bool
            Whether to add dispersion correction. Default is False.
    '''
    def __init__(self, mols, x, T, sigmas, parameters):
        self.mols = mols
        self.x = np.asarray(x)
        self.T = T
        self.sigmas = sigmas
        self.parameters = parameters
        self.split_sigma = sigmas[0].split_sigma
        self.dispersion = False
        self.sanity_check()

    def sanity_check(self):
        grid = self.sigmas[0].sigma_grid
        for i in range(1,len(self.sigmas)):
            if not np.allclose(grid, self.sigmas[i].sigma_grid):
                raise RuntimeError("inconsitent sigma grids in sigma files")

    def kernel(self):
        lngamma = lngamma_c(self.mols, self.x, self.parameters.parameters)
        if not self.split_sigma:
            lngamma += lngamma_r(self.mols, self.sigmas, self.x, self.T, self.parameters.parameters)
        else:
            lngamma += lngamma_r3(self.mols, self.sigmas, self.x, self.T, self.parameters.parameters)
        if self.dispersion:
            lngamma += lngamma_dsp(self.mols, self.sigmas, self.x, self.parameters.parameters)
        return lngamma


if __name__ == "__main__":
    from pycosmosac.param import parameters, data
    from pycosmosac.cosmo import cosmo
    from pycosmosac.sigma import sigma

    myparam = parameters.Parameters()

    mol1 = cosmo.Cosmo().load("./test/butane.cosmo")
    sigma1 = sigma.Sigma(mol1, myparam)
    sigma1.write_sigma_file = False
    sigma1.kernel()

    mol2 = cosmo.Cosmo().load("./test/h2o.cosmo")
    sigma2 = sigma.Sigma(mol2, myparam)
    sigma2.write_sigma_file = False
    sigma2.kernel()

    x = [0., 1.]
    T = 298.15
    myac = AC([mol1,mol2], x, T, [sigma1,sigma2], myparam)
    print(myac.kernel() - np.asarray([10.05122252, 0.0]))

    myparam = parameters.Parameters(data.Hsieh_2010)
    sigma1 = sigma.Sigma(mol1, myparam)
    sigma1.write_sigma_file = False
    sigma1.split_sigma = True
    sigma1.kernel()

    sigma2 = sigma.Sigma(mol2, myparam)
    sigma2.write_sigma_file = False
    sigma2.split_sigma = True
    sigma2.kernel()
    myac = AC([mol1,mol2], x, T, [sigma1,sigma2], myparam)
    print(myac.kernel() - np.asarray([8.81631154, 0.0]))

    myac.dispersion = True
    print(myac.kernel() - np.asarray([9.55918891, 0.]))
