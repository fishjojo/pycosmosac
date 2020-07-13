from .cavity import Cavity

class Mole():
    '''
    Class for molecular information
    '''
    def __init__(self):
        #{"atom" : [], "xyz" : ndarray(natom, 3)}
        self.geometry = {}
        self.cavity = Cavity()

    def build(self, geometry=None, cavity=None):
        if geometry is not None: self.geometry = geometry
        if cavity is not None: self.cavity = cavity
        return self
