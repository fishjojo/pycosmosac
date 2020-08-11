import os
import warnings
import numpy as np
from scipy.spatial import distance
from pycosmosac.molecule.cavity import Cavity
from pycosmosac.param import data
from pycosmosac.utils import elements

BOND_SCALING = 1.2

def get_connectivity(mol, geometry=None):
    #TODO improve accuracy
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
            warnings.warn("atom (%s, %s) has no bonds." % (i+1, atom_i))
    return connectivity

def _dfs(connectivity, iatm, color, traversalOrder, res, parent=None):
    '''
    depth-first search
    '''
    color[iatm] = 1
    traversalOrder.append(iatm)
    for jatm in connectivity[iatm]:
        if color[jatm] == 0:
            if len(connectivity[jatm]) < 2:
                color[jatm] = 2
            else:
                _dfs(connectivity, jatm, color, traversalOrder, res, parent=iatm)
        elif color[jatm] == 1:
            if parent and parent != jatm:
                cycle = []
                lastatm_index = traversalOrder.index(iatm)
                for index in range(lastatm_index, -1, -1):
                    katm = traversalOrder[index]
                    if katm == jatm:
                        break
                    else:
                        cycle.append(katm)
                cycle.append(jatm)
                res.append(cycle)
    color[iatm] = 2
    traversalOrder.pop()

def find_rings(mol, connectivity=None):
    if connectivity is None: connectivity = mol.connectivity
    natom = mol.natom
    color = np.zeros((natom), dtype=int)
    res = []
    for i in range(natom):
        if color[i] > 0:
            continue
        if len(connectivity[i]) < 2:
            color[i] = 2
            continue
        traversalOrder = []
        _dfs(connectivity, i, color, traversalOrder, res)
    return res

def find_ring_atoms(mol, connectivity=None):
    if connectivity is None: connectivity = mol.connectivity
    res = find_rings(mol, connectivity)
    return list(set().union(*res))

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
        elif atom_i == "Br":
            natom += 1
            disp_tot += data.disp["Br"]
        elif atom_i == "I":
            natom += 1
            disp_tot += data.disp["I"]
        elif atom_i == "P":
            natom += 1
            disp_tot += data.disp["P"]
        elif atom_i == "S":
            natom += 1
            disp_tot += data.disp["S"]
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

def fromstring(string, format='xyz'):
    format = format.lower()
    if format == 'xyz':
        dat = string.splitlines()
        natm = int(dat[0])
        return '\n'.join(dat[2:natm+2])
    elif format == 'raw':
        return string
    else:
        raise NotImplementedError

def fromfile(filename, format=None):
    if format is None:  # Guess format based on filename
        format = os.path.splitext(filename)[1][1:].lower()
        if format not in ('xyz', 'zmat', 'sdf', 'mol2'):
            format = 'raw'
    with open(filename, 'r') as f:
        return fromstring(f.read(), format)

def read_geometry(xyz):
    if os.path.isfile(xyz):
        try:
            xyz_raw = fromfile(xyz)
            return raw_to_geometry(xyz_raw)
        except:
            raise ValueError('Failed to parse geometry file  %s' % xyz)
    else:
        return raw_to_geometry(xyz)

def raw_to_geometry(xyz):
    geometry = {}
    geometry["atom"] = []
    geometry["xyz"] = []

    def str2atm(line):
        dat = line.split()
        assert(len(dat) == 4)
        geometry["atom"].append(dat[0])
        geometry["xyz"].append([float(x) for x in dat[1:4]])

    if isinstance(xyz, str):
        xyz = str(xyz.replace(';','\n').replace(',',' ').replace('\t',' '))
        fmt_atoms = []
        for dat in xyz.split('\n'):
            dat = dat.strip()
            if dat and dat[0] != '#':
                fmt_atoms.append(dat)
        for line in fmt_atoms:
            str2atm(line)
        geometry["xyz"] = np.asarray(geometry["xyz"])
    else:
        raise NotImplementedError
    return geometry

def geometry_to_xyz(geometry, name="unknown"):
    symb = geometry["atom"]
    coord = geometry["xyz"]
    natom = len(symb)
    xyz = str(natom) + "\n"
    xyz += name + "\n"
    for i in range(natom):
        xyz += symb[i] + "    "
        xyz += str(coord[i, 0]) + "    "
        xyz += str(coord[i, 1]) + "    "
        xyz += str(coord[i, 2]) + "    "
        xyz += "\n"
    return xyz.strip()

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

    @property
    def natom(self):
        if self.geometry is None:
            raise RuntimeError("molecule not initialized")
        return len(self.geometry["atom"])

    def build(self, geometry=None, cavity=None):
        if geometry is not None: 
            if isinstance(geometry, str):
                self.geometry = read_geometry(geometry)
            elif isinstance(geometry, dict):
                self.geometry = geometry
            else:
                raise ValueError("unsupported geometry input")
        if cavity is not None: self.cavity = cavity
        self.connectivity = self.get_connectivity()
        self.hb_class = self.classify_hydrogen_bonds()
        return self

    get_connectivity = get_connectivity
    classify_hydrogen_bonds = classify_hydrogen_bonds
    get_dispersion_type = get_dispersion_type
    find_rings = find_rings
    find_ring_atoms = find_ring_atoms

if __name__ == "__main__":
    from pycosmosac.utils.misc import fingerprint
    geometry = {}
    geometry["atom"] = ['O', 'H', 'H']
    geometry["xyz"] =  np.asarray([[ 0., 0., -0.405655705],
                                  [ 0.770106178, 0., 0.202827852],
                                  [-0.770106178, 0., 0.202827852]])
    mol = Mole().build(geometry)
    print(mol.connectivity == [[1, 2], [0], [0]])
    print(mol.hb_class == ['OH','OH','OH'])
    print(mol.get_dispersion_type()[0] - 70.75953333333332)
    print(mol.get_dispersion_type()[1] == "H2O")

    xyz = '''
        N         -2.86237        0.53549       -0.00680
        C         -1.59157        1.12789       -0.00460
        C         -0.65647        0.06499       -0.00640
        N         -1.36327       -1.15231       -0.01640
        C         -2.67117       -0.86471       -0.01510
        C          0.72143        0.40079       -0.00170
        N          1.06103        1.72159       -0.02520
        C          0.08733        2.69019       -0.01940
        N         -1.23787        2.45869       -0.00720
        H          0.42943        3.73339       -0.02780
        N          1.75063       -0.53271        0.12520
        H         -3.73057        1.00639       -0.00350
        H         -3.47277       -1.60891       -0.01910
        H          1.51683       -1.44251       -0.19520
        H          2.65133       -0.22311       -0.15800
    '''
    mol = Mole().build(xyz)
    print(mol.geometry["atom"] == ['N', 'C', 'C', 'N', 'C', 'C', 'N', 'C', 'N', 'H', 'N', 'H', 'H', 'H', 'H'])
    print(fingerprint(mol.geometry["xyz"]) - -7.705571225872962)
    print(mol.connectivity == [[1, 4, 11], [0, 2, 8], [1, 3, 5], [2, 4], [0, 3, 12], [2, 6, 10], [5, 7], [6, 8, 9], [1, 7], [7], [5, 13, 14], [0], [4], [10], [10]])
    print(mol.find_ring_atoms() == [0, 1, 2, 3, 4, 5, 6, 7, 8])
    print(mol.find_rings() == [[4, 3, 2, 1, 0], [8, 7, 6, 5, 2, 1]])

    xyz = '''
        C          0.78526        0.09180       -0.07290
        C          0.46366       -1.44010        0.02370
        C         -0.58284       -1.14900        1.15480
        C         -0.26134        0.38280        1.05820
        C         -0.33764        0.25930       -1.15480
        C         -1.38414        0.55030       -0.02370
        C         -1.70564       -0.98160        0.07290
        C         -0.65914       -1.27260       -1.05820
        H          1.78286        0.52170       -0.13130
        H          1.20356       -2.23730        0.04260
        H         -0.68114       -1.71320        2.07970
        H         -0.10194        1.04580        1.90590
        H         -0.23944        0.82320       -2.07990
        H         -2.12424        1.34730       -0.04270
        H         -2.70324       -1.41140        0.13130
        H         -0.81854       -1.93550       -1.90580
    '''
    mol = Mole().build(xyz)
    print(mol.find_ring_atoms() == [0, 1, 2, 3, 4, 5, 6, 7])
