import unittest
from pycosmosac.utils import antoine

class KnownValues(unittest.TestCase):
    def test_antoine_nist(self):
        data = antoine.get_antoine_nist("7732-18-5")
        P = antoine.antoine_to_vapor(data, 300)
        self.assertAlmostEqual(P, 3554.6837736037555, 9)

    def test_vapor(self):
        data = [[[250.04, 328.57], [350.14, 466.73], [212.4, 293.02]], 
                [4.022, 4.46988, 4.13377], 
                [1062.64, 1354.913, 1102.878], 
                [-44.93, -5.537, -40.46]]
        P = antoine.antoine_to_vapor(data, 298.15)
        self.assertAlmostEqual(P, 66910.00086027385, 9)

if __name__ == '__main__':
    unittest.main()

