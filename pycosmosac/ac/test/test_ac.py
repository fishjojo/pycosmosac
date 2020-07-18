import os
import unittest
import numpy as np
from pycosmosac.param import parameters, data
from pycosmosac.cosmo import cosmo
from pycosmosac.sigma import sigma
from pycosmosac.ac import ac

here = os.path.abspath(os.path.dirname(__file__))
mol1 = cosmo.Cosmo().load(here+"/butane.cosmo")
mol2 = cosmo.Cosmo().load(here+"/h2o.cosmo")
x = [0., 1.]
T = 298.15

class KnownValues(unittest.TestCase):
    def test_infinity_dilution(self):
        myparam = parameters.Parameters()
        sigma1 = sigma.Sigma(mol1, myparam)
        sigma1.write_sigma_file = False
        sigma1.kernel()
        sigma2 = sigma.Sigma(mol2, myparam)
        sigma2.write_sigma_file = False
        sigma2.kernel()

        myac = ac.AC([mol1,mol2], x, T, [sigma1,sigma2], myparam)
        self.assertTrue(np.allclose(myac.kernel(), np.asarray([10.05122252, 0.0])))

    def test_infinity_dilution3(self):
        myparam = parameters.Parameters(data.Hsieh_2010)
        sigma1 = sigma.Sigma(mol1, myparam)
        sigma1.write_sigma_file = False
        sigma1.split_sigma = True
        sigma1.kernel()
        sigma2 = sigma.Sigma(mol2, myparam)
        sigma2.write_sigma_file = False
        sigma2.split_sigma = True
        sigma2.kernel()
        myac = ac.AC([mol1,mol2], x, T, [sigma1,sigma2], myparam)
        self.assertTrue(np.allclose(myac.kernel(), np.asarray([8.81631154, 0.0])))

if __name__ == '__main__':
    unittest.main()
