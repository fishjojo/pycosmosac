import unittest
import numpy as np
from pycosmosac.molecule import mole

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

if __name__ == '__main__':
    unittest.main()
