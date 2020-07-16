import os
import unittest
import numpy as np
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
        pA = []
        with open("/tmp/tmp.sigma", 'r') as f:
            contents = f.read().splitlines()
        for i, line in enumerate(contents[3:]):
            pA.append(float(line.split()[1]))
        self.assertAlmostEqual(misc.fp(np.asarray(pA)), -1.917622937152116, 8)

        mysigma.write_sigma_file = True
        mysigma.sigma_filename = "/tmp/tmp1.sigma"
        mysigma.kernel()
        with open("/tmp/tmp1.sigma", 'r') as f:
            contents1 = f.read().splitlines()
        self.assertEqual(contents1, contents)

if __name__ == '__main__':
    unittest.main()
