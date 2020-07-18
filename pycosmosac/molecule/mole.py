import warnings
import numpy as np
from scipy.spatial import distance
from pycosmosac.molecule.cavity import Cavity
from pycosmosac.param import data
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

def get_dispersion_type(mol, geometry=None, connectivity=None):
    if geometry is None: geometry = mol.geometry
    if connectivity is None: connectivity = mol.connectivity
    atoms = geometry["atom"]

    if len(atoms) == 3 and atoms.count("O") == 1 and atoms.count("H") == 2:
        disp_tot = (data.disp["H(H2O)"] * 2 + data.disp["-O-"]) / 3.0
        return disp_tot, "H2O"

    disp_type = "NHB"
    disp_tot = 0.0
    natom = 0
    nCOOH = 0
    for i, atom_i in enumerate(atoms):
        n = len(connectivity[i])
        if atom_i == "C":
            natom += 1
            if n == 4:
                disp_tot += data.disp["C(sp3)"] 
            elif n == 3:
                disp_tot += data.disp["C(sp2)"]
                atom_js = []
                js = []
                for j in connectivity[i]:
                    atom_js.append(atoms[j])
                    js.append(j)
                if atom_js.count("O") == 2:
                    for j, atom_j in zip(js,atom_js):
                        if atom_j != "O":
                            continue
                        if len(connectivity[j]) == 2:
                            for k in connectivity[j]:
                                if atoms[k] == "H":
                                    nCOOH += 1
                                    disp_tot += data.disp["H(COOH)"]
                                    disp_type = "COOH"
            elif n == 2:
                disp_tot += data.disp["C(sp)"]
        elif atom_i == "N":
            natom += 1
            if n == 3:
                disp_tot += data.disp["N(sp3)"]
            elif n == 2:
                disp_tot += data.disp["N(sp2)"]
            elif n == 1:
                disp_tot += data.disp["N(sp)"]
        elif atom_i == "O":
            natom += 1
            if n == 2:
                disp_tot += data.disp["-O-"]
            elif n == 1:
                disp_tot += data.disp["=O"]
        elif atom_i == "F":
            natom += 1
            disp_tot += data.disp["F"]
        elif atom_i == "Cl":
            natom += 1
            disp_tot += data.disp["Cl"]
        elif atom_i == "H":
            j = connectivity[i][0]
            atom_j = atoms[j]
            if atom_j == "O":
                natom += 1
                disp_tot += data.disp["H(OH)"]
            elif atom_j == "N":
                natom += 1
                disp_tot += data.disp["H(NH)"]
        else:
            warnings.warn("dispersion parameter not available for %s" % atom_i)

    disp_tot -= nCOOH * data.disp["H(OH)"]
    disp_tot /= natom

    return disp_tot, disp_type


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
    get_dispersion_type = get_dispersion_type

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
    print(mol.get_dispersion_type()[0] - 70.75953333333332)
    print(mol.get_dispersion_type()[1] == "H2O")
