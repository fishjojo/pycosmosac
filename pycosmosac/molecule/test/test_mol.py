import os
import unittest
import numpy as np
from pycosmosac.molecule import mole
from pycosmosac.cosmo import cosmo

here = os.path.abspath(os.path.dirname(__file__))
filename = here + "/150-13-0.cosmo"

geometry = {}
geometry["atom"] = ['O', 'H', 'H']
geometry["xyz"] =  np.asarray([[ 0., 0., -0.405655705],
                               [ 0.770106178, 0., 0.202827852],
                               [-0.770106178, 0., 0.202827852]])

class KnownValues(unittest.TestCase):
    def test_build(self):
        mol = mole.Mole()
        mol.build(geometry)
        self.assertEqual(mol.connectivity, [[1, 2], [0], [0]])
        self.assertEqual(mol.hb_class, ['OH','OH','OH'])

    def test_dsp(self):
        mycosmo = cosmo.Cosmo()
        mycosmo.load(filename)
        mol = mycosmo.mol
        disp_eps, disp_type = mol.get_dispersion_type()
        self.assertAlmostEqual(disp_eps, 97.15234615384617, 9)
        self.assertEqual(disp_type, "COOH")

if __name__ == '__main__':
    unittest.main()
