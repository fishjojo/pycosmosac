import numpy as np

class Cavity():
    '''
    Class for molecular cavity
    '''
    def __init__(self):
        self.area = None
        self.volume = None
        #{"xyz":ndarray(n,3), "charge":ndarray(n,), "area":ndarray(n,), "charge/area":ndarray(n,)}
        self.segments = {}
        self.atom_map = np.empty((0), dtype=int)
