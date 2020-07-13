import re
from io import StringIO
import numpy as np
import pandas
from pycosmosac.molecule import mole, cavity

def get_molecule(data):
    sdata = re.search(r"!DATE[a-zA-Z0-9:\s]+\n(.+)end\s*", data, re.DOTALL).group(1)
    df = pandas.read_csv(StringIO(sdata), names=['atomidentifier','x / A','y / A','z / A','?1','?2','?3','atom','?4'],sep=r'\s+',engine = 'python')

    geometry = {}
    x = np.asarray(df['x / A'].tolist(), dtype=float)
    y = np.asarray(df['y / A'].tolist(), dtype=float)
    z = np.asarray(df['z / A'].tolist(), dtype=float)
    geometry["xyz"] = np.column_stack((x,y,z))
    geometry["atom"] = df['atom'].tolist()
    return geometry

def get_cavity(data):
    cav = cavity.Cavity()
    cav.area = float(re.search(r"Total surface area of cavity \(A\*\*2\)     =(.+)\n", data).group(1).strip())
    cav.volume = float(re.search(r"Total volume of cavity \(A\*\*3\)           =(.+)\n", data).group(1).strip())

    sdata = re.search(r"\(X, Y, Z\)[\sa-zA-Z0-9\[\]\^\./]+\n(.+)(\n\n|$)", data, re.DOTALL).group(1).rstrip()
    df = pandas.read_csv(StringIO(sdata), names=['n','atom','x / a.u.','y / a.u.','z / a.u.','charge / e','area / A^2','charge/area / e/A^2','potential'],sep=r'\s+',engine= 'python')
   
    cav.atom_map = np.asarray(df['atom'].tolist(), dtype=int)
    x = np.asarray(df['x / a.u.'].tolist(), dtype=float)
    y = np.asarray(df['y / a.u.'].tolist(), dtype=float)
    z = np.asarray(df['z / a.u.'].tolist(), dtype=float)
    cav.segments["xyz"] = np.column_stack((x,y,z))
    cav.segments["charge"] = np.asarray(df['charge / e'].tolist(), dtype=float)
    cav.segments["area"] = np.asarray(df['area / A^2'].tolist(), dtype=float)
    cav.segments["charge/area"] = np.asarray(df['charge/area / e/A^2'].tolist(), dtype=float)
    return cav

def load(name):
    try:
        with open(name, 'r') as f:
            data = f.read()
            geometry = get_molecule(data)
            cavity = get_cavity(data)
            mol = mole.Mole()
            mol.build(geometry = geometry, cavity = cavity)
            return mol
    except:
        raise RuntimeError("failed reading cosmo file %s." % name)


if __name__ == "__main__":
    mol = load("h2o.cosmo")
    print(mol.geometry)
    #print(mol.cavity.segments)
