
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
    import numpy as np

    mycosmo = Cosmo()
    mycosmo.load("./test/h2o.cosmo")
    cavity = mycosmo.mol.cavity
    total_area = cavity.area
    segments = cavity.segments
    charge = segments["charge"]
    area = segments["area"]
    print(np.sum(area) - total_area)

