import numpy as np
from scipy.spatial import distance
from pycosmosac.molecule.cavity import Cavity
from pycosmosac.utils import elements

BOND_SCALING = 1.2

def get_connectivity(mol, geometry=None):
    if geometry is None: 
        geometry = mol.geometry
    if not geometry:
        raise RuntimeError("molecule not initialized.")
    atoms = geometry["atom"]
    xyz = geometry["xyz"]
    d = distance.cdist(xyz, xyz)

    natom = len(atoms)
    connectivity = [[] for _ in range(natom)]
    for i, atom_i in enumerate(atoms):
        for j, atom_j in enumerate(atoms):
            if i==j:
                continue
            l = BOND_SCALING * elements.covalent_bond(atom_i, atom_j)
            if d[i,j] <= l:
                connectivity[i].append(j)
        if not connectivity[i]:
            raise RuntimeError("atom (%s, %s) has no bonds." % (i+1, atom_i))
    return connectivity

def classify_hydrogen_bonds(mol, geometry=None, connectivity=None):
    if geometry is None: geometry = mol.geometry
    if connectivity is None: connectivity = mol.connectivity
    if not geometry or not connectivity:
        raise RuntimeError("molecule not initialized.")

    atoms = geometry["atom"]
    hb_class = []
    for i, atom_i in enumerate(atoms):
        if atom_i in ['N', 'F']:
            hb_class.append("OT")
        elif atom_i in ['O', 'H']:
            bond_type = 'NHB'
            for j in connectivity[i]:
                atom_j = atoms[j]
                atom_ij = atom_i + atom_j
                if  atom_ij in ['OH', 'HO']:
                    bond_type = 'OH'
                    break
                if atom_i == 'O':
                    bond_type = 'OT'
                    break
                if atom_i == 'H' and atom_j in ['N', 'F']:
                    bond_type = 'OT'
                    break
            hb_class.append(bond_type)
        else:
            hb_class.append('NHB')
    return hb_class

class Mole():
    '''
    Class for molecular information
    Attributes:
        geometry : dict
            Geometry information.
                "xyz" : ndarray
                "atom" : list
        cavity : Cavity
            Cavity information
        connectivity : list
            Connectivity information.
        hb_class : list
            Hydrogen bond classification.
    '''
    def __init__(self):
        #{"atom" : [], "xyz" : ndarray(natom, 3)}
        self.geometry = None
        self.cavity = None
        self.connectivity = None
        self.hb_class = None

    def build(self, geometry=None, cavity=None):
        if geometry is not None: self.geometry = geometry
        if cavity is not None: self.cavity = cavity
        self.connectivity = self.get_connectivity()
        self.hb_class = self.classify_hydrogen_bonds()
        return self

    get_connectivity = get_connectivity
    classify_hydrogen_bonds = classify_hydrogen_bonds

if __name__ == "__main__":
    mol = Mole()
    geometry = {}
    geometry["atom"] = ['O', 'H', 'H']
    geometry["xyz"] =  np.asarray([[ 0., 0., -0.405655705],
                                  [ 0.770106178, 0., 0.202827852],
                                  [-0.770106178, 0., 0.202827852]])
    mol.build(geometry)
    print(mol.connectivity == [[1, 2], [0], [0]])
    print(mol.hb_class == ['OH','OH','OH'])
