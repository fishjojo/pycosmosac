import warnings
import numpy as np
from pycosmosac.utils import constants as const
from pycosmosac.ac import ac

SOLVE_GAMMA_THRESH = 1e-3
SOLVE_GAMMA_MAXITER = 200

def mu_c(mols, x, parameters):
    '''
    Combinatorial contribution based on COSMOtherm release 15.01
    '''
    r0 = parameters["r_0"]
    q0 = parameters["q_0"]
    lambda0 = parameters["lambda_0"]
    lambda1 = parameters["lambda_1"]
    lambda2 = parameters["lambda_2"]

    A = []
    V = []
    for mol in mols:
        A.append(mol.cavity.area)
        V.append(mol.cavity.volume)
    A = np.asarray(A)
    V = np.asarray(V)

    r = V / r0
    q = A / q0

    xq = np.sum(x * q)
    xr = np.sum(x * r)
    r_xr = r / xr
    q_xq = q / xq

    RT = const.R * T
    muc = RT * (lambda0 * np.log(r) + lambda1 * (1.0 - r_xr + np.log(r_xr)) + lambda2 * (1.0 - q_xq + np.log(q_xq)))
    return muc


def solve_sigma_potential(W, ps, aeff, T, thresh=SOLVE_GAMMA_THRESH, maxiter=SOLVE_GAMMA_MAXITER):
    mus = np.zeros_like(ps)
    RT = const.R * T
    count = 0
    while (count < maxiter):
        mus_new = -RT / aeff * np.log(np.sum(np.exp(aeff * (mus - W) / RT) * ps, axis=1))
        mus = 0.5 * (mus + mus_new)
        diff = np.amax(np.absolute(mus_new - mus))
        if diff < thresh:
            break
        count += 1
    if count == maxiter:
        warnings.warn("solving sigma potential reached max iterations: %s steps with diff %.3e" % (maxiter,diff))
    return mus


def mu_s(mols, sigmas, x, T, parameters, thresh=SOLVE_GAMMA_THRESH, maxiter=SOLVE_GAMMA_MAXITER):
    aeff = parameters["a_eff"]
    A = []
    for mol in mols:
        A.append(mol.cavity.area)
    A = np.asarray(A)

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
    W = ac.calc_eint(grid, alpha, sigma_hb, chb)
    mus = solve_sigma_potential(W, ps, aeff, T, thresh, maxiter)

    result = []
    for pi in p:
        result.append(np.dot(pi, mus))
    result = np.asarray(result)
    return result


class RS(ac.AC):

    def __init__(self, mols, x, T, sigmas, parameters):
        ac.AC.__init__(self, mols, x, T, sigmas, parameters)
        self.comb = True

    def kernel(self):
        mu = mu_s(self.mols, self.sigmas, self.x, self.T, self.parameters.parameters)
        if self.comb:
            mu += mu_c(self.mols, self.x, self.parameters.parameters)
        return mu


if __name__ == "__main__":
    from pycosmosac.param import parameters, data
    from pycosmosac.cosmo import cosmo
    from pycosmosac.sigma import sigma

    myparam = parameters.Parameters(data.COSMORS)

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
    myac = RS([mol1,mol2], x, T, [sigma1,sigma2], myparam)
    print(myac.kernel())

