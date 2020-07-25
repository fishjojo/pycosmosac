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

ATOMIC_MASS = [
    0.,                 # GHOST
    1.008,              # H [1.00784, 1.00811]
    4.002602,           # He
    6.94,               # Li [6.938, 6.997]
    9.0121831,          # Be
    10.81,              # B [10.806, 10.821]
    12.011,             # C [12.0096, 12.0116]
    14.007,             # N [14.00643, 14.00728]
    15.999,             # O [15.99903, 15.99977]
    18.998403163,       # F
    20.1797,            # Ne
    22.98976928,        # Na
    24.305,             # Mg [24.304, 24.307]
    26.9815385,         # Al
    28.085,             # Si [28.084, 28.086]
    30.973761998,       # P
    32.06,              # S [32.059, 32.076]
    35.45,              # Cl [35.446, 35.457]
    39.948,             # Ar
    39.0983,            # K
    40.078,             # Ca
    44.955908,          # Sc
    47.867,             # Ti
    50.9415,            # V
    51.9961,            # Cr
    54.938044,          # Mn
    55.845,             # Fe
    58.933194,          # Co
    58.6934,            # Ni
    63.546,             # Cu
    65.38,              # Zn
    69.723,             # Ga
    72.630,             # Ge
    74.921595,          # As
    78.971,             # Se
    79.904,             # Br [79.901, 79.907]
    83.798,             # Kr
    85.4678,            # Rb
    87.62,              # Sr
    88.90584,           # Y
    91.224,             # Zr
    92.90637,           # Nb
    95.95,              # Mo
    97.90721,           # 98Tc
    101.07,             # Ru
    102.90550,          # Rh
    106.42,             # Pd
    107.8682,           # Ag
    112.414,            # Cd
    114.818,            # In
    118.710,            # Sn
    121.760,            # Sb
    127.60,             # Te
    126.90447,          # I
    131.293,            # Xe
    132.90545196,       # Cs
    137.327,            # Ba
    138.90547,          # La
    140.116,            # Ce
    140.90766,          # Pr
    144.242,            # Nd
    144.91276,          # 145Pm
    150.36,             # Sm
    151.964,            # Eu
    157.25,             # Gd
    158.92535,          # Tb
    162.500,            # Dy
    164.93033,          # Ho
    167.259,            # Er
    168.93422,          # Tm
    173.054,            # Yb
    174.9668,           # Lu
    178.49,             # Hf
    180.94788,          # Ta
    183.84,             # W
    186.207,            # Re
    190.23,             # Os
    192.217,            # Ir
    195.084,            # Pt
    196.966569,         # Au
    200.592,            # Hg
    204.38,             # Tl [204.382, 204.385]
    207.2,              # Pb
    208.98040,          # Bi
    208.98243,          # Po
    209.98715,          # At
    222.01758,          # Rn
    223.01974,          # Fr
    226.02541,          # Ra
    227.02775,          # Ac
    232.0377,           # Th
    231.03588,          # Pa
    238.02891,          # U
    237.04817,          # Np
    244.06421,          # Pu
    243.06138,          # Am
    247.07035,          # Cm
    247.07031,          # Bk
    251.07959,          # Cf
    252.0830,           # Es
    257.09511,          # Fm
    258.09843,          # Md
    259.1010,           # No
    262.110,            # Lr
    267.122,            # Rf
    268.126,            # Db
    271.134,            # Sg
    270.133,            # Bh
    269.1338,           # Hs
    278.156,            # Mt
    281.165,            # Ds
    281.166,            # Rg
    285.177,            # Cn
    286.182,            # Nh
    289.190,            # Fl
    289.194,            # Mc
    293.204,            # Lv
    293.208,            # Ts
    294.214,            # Og
]

def molecular_mass(mol):
    atoms = mol.geometry["atom"]
    mass = 0.0
    for a in atoms:
        ia = ELEMENTS.index(a)
        mass += ATOMIC_MASS[ia]
    return mass
