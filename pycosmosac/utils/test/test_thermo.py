import unittest
from pycosmosac.utils import thermo

class KnownValues(unittest.TestCase):
    def test_G_binary(self):
        lngamma = 3.0
        T = 300
        G = thermo.calc_G_binary(lngamma, T)
        self.assertAlmostEqual(G, 1.7884858072629, 9)

        density = [998.]
        molar_mass = [18.]
        vaporP = 1e4
        G = thermo.calc_G_binary(lngamma, T, vaporP=vaporP, density=density, molar_mass=molar_mass)
        self.assertAlmostEqual(G, -3.8956651094414827, 9)

        density = [800., 998.]
        molar_mass = [20., 18.]
        vaporP = 1e4
        P = 1e5
        G = thermo.calc_G_binary(lngamma, T, vaporP=vaporP, P=P, density=density, molar_mass=molar_mass)
        self.assertAlmostEqual(G, -3.8951273459414826, 9)

if __name__ == '__main__':
    unittest.main()
