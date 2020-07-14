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

def solve_lnGamma(W, GammaS, ps, T, thresh=1e-6, maxiter=500):
    #TODO improve convergence speed
    count = 0
    while (count < maxiter):
        GammaS_new = np.exp(-np.log(np.sum(np.exp(-W / (const.R*T)) * GammaS * ps, axis=1)))
        GammaS = 0.5 * (GammaS + GammaS_new)
        diff = np.amax(np.absolute(GammaS_new - GammaS))
        if diff < thresh:
            break
        count += 1
    return np.log(GammaS)

def lngamma_r(mols, sigmas, x, T, parameters, thresh=1e-6, maxiter=500):
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
    '''
    def __init__(self, mols, x, T, sigmas, parameters):
        self.mols = mols
        self.x = np.asarray(x)
        self.T = T
        self.sigmas = sigmas
        self.parameters = parameters
        self.sanity_check()

    def sanity_check(self):
        grid = self.sigmas[0].sigma_grid
        for i in range(1,len(self.sigmas)):
            if not np.allclose(grid, self.sigmas[i].sigma_grid):
                raise RuntimeError("inconsitent sigma grids in sigma files")

    def kernel(self):
        lngamma =  lngamma_c(self.mols, self.x, self.parameters.parameters)
        lngamma += lngamma_r(self.mols, self.sigmas, self.x, self.T, self.parameters.parameters)
        return lngamma


if __name__ == "__main__":
    from pycosmosac.param import parameters
    from pycosmosac.cosmo import cosmo
    from pycosmosac.sigma import sigma

    myparam = parameters.Parameters()

    mol1 = cosmo.Cosmo().load("butane.cosmo")
    sigma1 = sigma.Sigma(mol1, myparam)
    sigma1.write_sigma_file = False
    sigma1.kernel()

    mol2 = cosmo.Cosmo().load("h2o.cosmo")
    sigma2 = sigma.Sigma(mol2, myparam)
    sigma2.write_sigma_file = False
    sigma2.kernel()

    x = [0., 1.]
    T = 298.15
    myac = AC([mol1,mol2], x, T, [sigma1,sigma2], myparam)
    print(myac.kernel())
