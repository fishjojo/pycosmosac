import itertools

ELEMENTS = ['X', #Ghost
    'H' , 'He', 'Li', 'Be', 'B' , 'C' , 'N' , 'O' , 'F' , 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P' , 'S' , 'Cl', 'Ar', 'K' , 'Ca',
    'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y' , 'Zr',
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
    'Sb', 'Te', 'I' , 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    'Lu', 'Hf', 'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
    'Pa', 'U' , 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
    'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
    'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og',
]

ELEMENTS_UPPER = dict((x.upper(), x) for x in ELEMENTS)

def std_symb(symb_or_chg):
    if isinstance(symb_or_chg, str):
        symb = symb_or_chg.upper()
        if symb in ELEMENTS_UPPER:
            return ELEMENTS_UPPER[symb]
        else:
            raise ValueError("unrecognized symbol: %s." % symb_or_chg)
    else:
        return ELEMENTS[symb_or_chg]

# Covalent radii in Angstrom (https://doi.org/10.1039/b801115j)
COVALENT_RADIUS = {
    'H':0.31,
    'He':0.28,
    'Li':1.28,
    'Be':0.96,
    'B':0.84,
    'C':0.76, # sp3 hybridization; sp2: 0.73; sp: 0.69
    'N':0.71,
    'O':0.66,
    'F':0.57,
    'Ne':0.58,
    'Na':1.66,
    'Mg':1.41,
    'Al':1.21,
    'Si':1.11,
    'P':1.07,
    'S':1.05,
    'Cl':1.02,
    'Ar':1.06,
    'K':2.03,
    'Ca':1.76,
    'Sc':1.70,
    'Ti':1.60,
    'V':1.53,
    'Cr':1.39,
    'Mn':1.39, # low spin; high spin: 1.61
    'Fe':1.32, # low spin; high spin: 1.52
    'Co':1.26, # low spin; high spin: 1.50
    'Ni':1.24,
    'Cu':1.32,
    'Zn':1.22,
    'Ga':1.22,
    'Ge':1.20,
    'As':1.19,
    'Se':1.20,
    'Br':1.20,
    'Kr':1.16,
    'Rb':2.20,
    'Sr':1.95,
    'Y':1.90,
    'Zr':1.75,
    'Nb':1.64,
    'Mo':1.54,
    'Tc':1.47,
    'Ru':1.46,
    'Rh':1.42,
    'Pd':1.39,
    'Ag':1.45,
    'Cd':1.44,
    'In':1.42,
    'Sn':1.39,
    'Sb':1.39,
    'Te':1.38,
    'I':1.39,
    'Xe':1.40,
    'Cs':2.44,
    'Ba':2.15,
    'La':2.07,
    'Ce':2.04,
    'Pr':2.03,
    'Nd':2.01,
    'Pm':1.99,
    'Sm':1.98,
    'Eu':1.98,
    'Gd':1.96,
    'Tb':1.94,
    'Dy':1.92,
    'Ho':1.92,
    'Er':1.89,
    'Tm':1.90,
    'Yb':1.87,
    'Lu':1.87,
    'Hf':1.75,
    'Ta':1.70,
    'W':1.62,
    'Re':1.51,
    'Os':1.44,
    'Ir':1.41,
    'Pt':1.36,
    'Au':1.36,
    'Hg':1.32,
    'Tl':1.45,
    'Pb':1.46,
    'Bi':1.48,
    'Po':1.40,
    'At':1.50,
    'Rn':1.50,
    'Fr':2.60,
    'Ra':2.21,
    'Ac':2.15,
    'Th':2.06,
    'Pa':2.00,
    'U':1.96,
    'Np':1.90,
    'Pu':1.87,
    'Am':1.80,
    'Cm':1.69
}

COVALENT_BONDS= {
    (i, j): COVALENT_RADIUS[i] + COVALENT_RADIUS[j] for i, j in itertools.product(COVALENT_RADIUS.keys(), repeat=2)
}

def covalent_bond(atom_i, atom_j):
    return COVALENT_BONDS[(atom_i, atom_j)]
