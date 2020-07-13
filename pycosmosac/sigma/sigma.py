

class Sigma():
    '''
    Class for sigma profiles
    '''
    def __init__(self):




    def write_sigma_file(self, name):
        try:
            with open(name, 'w') as f:
        except:
            raise RuntimeError("failed to write sigma file.")
