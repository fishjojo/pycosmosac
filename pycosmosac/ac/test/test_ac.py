import os
import unittest
import numpy as np
from pycosmosac.param import parameters
from pycosmosac.cosmo import cosmo
from pycosmosac.sigma import sigma
from pycosmosac.ac import ac

here = os.path.abspath(os.path.dirname(__file__))

myparam = parameters.Parameters()
mol1 = cosmo.Cosmo().load(here+"/butane.cosmo")
sigma1 = sigma.Sigma(mol1, myparam)
sigma1.write_sigma_file = False
sigma1.kernel()
mol2 = cosmo.Cosmo().load(here+"/h2o.cosmo")
sigma2 = sigma.Sigma(mol2, myparam)
sigma2.write_sigma_file = False
sigma2.kernel()
T = 298.15

class KnownValues(unittest.TestCase):
    def test_infinity_dilution(self):
        x = [0., 1.]
        myac = ac.AC([mol1,mol2], x, T, [sigma1,sigma2], myparam)
        self.assertTrue(np.allclose(myac.kernel(), np.asarray([10.05122252, 0.0])) == True)

if __name__ == '__main__':
    unittest.main()
