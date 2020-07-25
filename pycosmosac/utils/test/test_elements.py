import unittest
import numpy as np
from pycosmosac.molecule import mole
from pycosmosac.utils import elements

mol = mole.Mole()
geometry = {}
geometry["atom"] = ['O', 'H', 'H']
geometry["xyz"] =  np.asarray([[ 0., 0., -0.405655705],
                               [ 0.770106178, 0., 0.202827852],
                               [-0.770106178, 0., 0.202827852]])
mol.build(geometry)

class KnownValues(unittest.TestCase):
    def test_molecular_mass(self):
        mass = elements.molecular_mass(mol)
        self.assertAlmostEqual(mass, 18.015, 9)

    def test_covalent_bond(self):
        O = mol.geometry["atom"][0]
        H1 = mol.geometry["atom"][1]
        bond_length = elements.covalent_bond(O, H1)
        self.assertAlmostEqual(bond_length, 0.97, 9)

    def test_std_symbol(self):
        atoms = mol.geometry["atom"]
        symbs = []
        for atom in atoms:
            symb = elements.std_symb(atom.lower())
            symbs.append(symb)
        self.assertEqual(symbs, ["O", "H", "H"])
            
if __name__ == '__main__':
    unittest.main()

