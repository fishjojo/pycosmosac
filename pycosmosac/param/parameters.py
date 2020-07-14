import json
from pycosmosac.param import data

KEYS = ["a_eff", "f_decay", "q_0", "r_0", "z"]

class Parameters():
    '''
    Parameters for COSMO-SAC models.
    Attributes:
        parameters : dict
            user input paramters
    '''

    def __init__(self, parameters=data.Saidi_2002):
        self.parameters = None
        self.load(parameters)
        self.sanity_check()

    def load(self, parameters):
        '''
        Load parameters from a dict or from a JSON file
        '''
        if isinstance(parameters, dict):
            self.parameters = parameters
        elif isinstance(parameters, str):
            try:
                with open(parameters, 'r') as f:
                    self.parameters = json.load(f)
            except:
                raise ValueError("JSON format expected.")
        else:
            raise TypeError("str or dict expected, while the input is %s." % type(parameters))

    def sanity_check(self):
        for key in KEYS:
            if not key in self.parameters:
                raise RuntimeError("parameter set has to include %s." % key)

    def dump_to_string(self):
        '''
        Dump parameters to string
        '''
        s = json.dumps(self.parameters)
        return s

    def dump_to_json(self, name):
        '''
        Dump parameters to JSON file
        '''
        try:
            with open(name,'w') as f:
                json.dump(self.parameters, f)
        except:
            raise RuntimeError("dumping parameters to JSON file failed.")

if __name__ == "__main__":
    myparam = Parameters()
    print(myparam.dump_to_string())
    myparam.dump_to_json('parameters.json')
