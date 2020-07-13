import json

class Parameters():
    '''
    Parameters for COSMO-SAC models.
    Attributes:
        parameters : dict
            user input paramters
    '''

    def __init__(self, parameters="Mullins_2006.json"):
        self.parameters = None
        self.load(parameters)

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
