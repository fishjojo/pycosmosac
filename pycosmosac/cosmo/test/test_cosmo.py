import os
import unittest
import numpy as np
from pycosmosac.cosmo import cosmo
from pycosmosac.utils import misc

here = os.path.abspath(os.path.dirname(__file__))
filename = here + "/h2o.cosmo"
mycosmo = cosmo.Cosmo()

class KnownValues(unittest.TestCase):
    def test_load(self):
        mycosmo.load(filename)
        mol = mycosmo.mol
        xyz = mol.geometry["xyz"]
        atom= mol.geometry["atom"]
        self.assertAlmostEqual(misc.fp(xyz), -1.3049969366762706, 8)
        self.assertEqual(atom, ['O', 'H', 'H'])
        cavity = mol.cavity
        total_area = cavity.area
        self.assertAlmostEqual(total_area, 43.223568, 8)
        total_volume = cavity.volume
        self.assertAlmostEqual(total_volume, 25.661981, 8)
        segments = cavity.segments
        charge = segments["charge"]
        area = segments["area"]
        self.assertAlmostEqual(np.sum(area), total_area, 6)
        sigma = charge / area
        sigma0 = segments["sigma"]
        self.assertTrue(np.allclose(sigma, sigma0) == True)
        xyz_seg = segments["xyz"]
        self.assertAlmostEqual(misc.fp(xyz_seg), 14.347472646595433, 8)

if __name__ == '__main__':
    unittest.main()
