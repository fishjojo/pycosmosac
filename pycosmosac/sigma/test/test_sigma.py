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

    def test_dump(self):
        mysigma.dump_to_file("/tmp/tmp.sigma")
        with open("/tmp/tmp.sigma", 'r') as f:
            contents = f.read().splitlines()
            sigma, pA = contents[12].split()
            self.assertAlmostEqual(float(sigma), -0.0160, 9)
            self.assertAlmostEqual(float(pA), 1.63558853113542, 9)

if __name__ == '__main__':
    unittest.main()
