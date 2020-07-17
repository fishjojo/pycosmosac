import os
import unittest
import numpy as np
from pycosmosac.param import parameters, data
from pycosmosac.cosmo import cosmo
from pycosmosac.sigma import sigma
from pycosmosac.utils import misc

here = os.path.abspath(os.path.dirname(__file__))

mycosmo = cosmo.Cosmo()
mol = mycosmo.load(here+"/h2o.cosmo")

param1 = parameters.Parameters()
sigma1 = sigma.Sigma(mol, param1)
sigma1.write_sigma_file = False
sigma1.kernel()

param3 = parameters.Parameters(data.Hsieh_2010)
sigma3 = sigma.Sigma(mol, param3)
sigma3.split_sigma = True
sigma3.write_sigma_file = False
sigma3.kernel()

class KnownValues(unittest.TestCase):
    def test_sigma(self):
        self.assertAlmostEqual(misc.fp(sigma1.pA), -1.917622937152116, 8)

    def test_sigma3(self):
        self.assertTrue(np.allclose(sigma3.pA, sigma3.pA_nhb + sigma3.pA_oh + sigma3.pA_ot))
        self.assertAlmostEqual(misc.fp(sigma3.pA_nhb), 0.049856600131970685, 8)
        self.assertAlmostEqual(misc.fp(sigma3.pA_oh), -0.6702525192899301, 8)
        self.assertAlmostEqual(misc.fp(sigma3.pA_ot), 0.0, 8)

    def test_dump(self):
        sigma1.dump_to_file("/tmp/tmp.sigma")
        pA = []
        with open("/tmp/tmp.sigma", 'r') as f:
            contents = f.read().splitlines()
        for i, line in enumerate(contents[3:]):
            pA.append(float(line.split()[1]))
        self.assertAlmostEqual(misc.fp(np.asarray(pA)), -1.917622937152116, 8)

        sigma1.write_sigma_file = True
        sigma1.sigma_filename = "/tmp/tmp1.sigma"
        sigma1.kernel()
        with open("/tmp/tmp1.sigma", 'r') as f:
            contents1 = f.read().splitlines()
        self.assertEqual(contents1, contents)

    def test_dump3(self):
        sigma3.dump_to_file("/tmp/tmp.sigma")
        pA = []
        with open("/tmp/tmp.sigma", 'r') as f:
            contents = f.read().splitlines()
        for i, line in enumerate(contents[3:]):
            pA.append(float(line.split()[1]))
        self.assertAlmostEqual(misc.fp(np.asarray(pA)), 0.8390277933438104, 8)

        sigma3.write_sigma_file = True
        sigma3.sigma_filename = "/tmp/tmp1.sigma"
        sigma3.kernel()
        with open("/tmp/tmp1.sigma", 'r') as f:
            contents1 = f.read().splitlines()
        self.assertEqual(contents1, contents)


if __name__ == '__main__':
    unittest.main()
