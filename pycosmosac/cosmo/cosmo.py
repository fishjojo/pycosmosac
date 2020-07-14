
class Cosmo():
    '''
    Class for reading cosmo files
    Attributes:
        mol : Mole
            Molecular information.
    '''
    def __init__(self):
        self.mol = None

    def load(self, filename, fn=None):
        if fn is None or not callable(fn):
            from pycosmosac.cosmo.read_cosmo_bagel import load as load_bagel
            fn = load_bagel
        self.mol = fn(filename)
        return self.mol

if __name__ == "__main__":
    mycosmo = Cosmo()
    mycosmo.load("h2o.cosmo")
    print(mycosmo.mol.geometry)
