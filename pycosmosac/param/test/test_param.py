import os
import unittest
import json
from pycosmosac.param import parameters, data

here = os.path.abspath(os.path.dirname(__file__))

class KnownValues(unittest.TestCase):
    def test_load(self):
        myparam = parameters.Parameters() #default
        self.assertEqual(myparam.parameters["a_eff"], 6.22, 9)

        myparam = parameters.Parameters(parameters=data.Mullins_2006)
        self.assertEqual(myparam.parameters["a_eff"], 7.5, 9)

        myparam = parameters.Parameters(parameters=here+"/Hsieh_2010.json")
        self.assertEqual(myparam.parameters["a_eff"], 7.25, 9)

    def test_dump(self):
        myparam = parameters.Parameters()
        s = myparam.dump_to_string()
        param = json.loads(s)
        self.assertEqual(param["a_eff"], 6.22, 9)

        myparam.dump_to_json("/tmp/tmp.json")
        with open("/tmp/tmp.json", 'r') as f:
            param = json.load(f)
            self.assertEqual(param["a_eff"], 6.22, 9)

if __name__ == '__main__':
    unittest.main()

