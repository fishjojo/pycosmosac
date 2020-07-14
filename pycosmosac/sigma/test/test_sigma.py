import os
import unittest
from pycosmosac.param import parameters
from pycosmosac.cosmo import cosmo
from pycosmosac.sigma import sigma
from pycosmosac.utils import misc

here = os.path.abspath(os.path.dirname(__file__))

myparam = parameters.Parameters()
mycosmo = cosmo.Cosmo()
mol = mycosmo.load(here+"/h2o.cosmo")
mysigma = sigma.Sigma(mol, myparam)
mysigma.write_sigma_file = False
mysigma.kernel()

class KnownValues(unittest.TestCase):
    def test_sigma(self):
        self.assertAlmostEqual(misc.fp(mysigma.pA), -1.917622937152116, 8)

if __name__ == '__main__':
    unittest.main()
