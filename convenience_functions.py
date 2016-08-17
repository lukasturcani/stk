import rdkit.Chem as chem
import rdkit.Chem.AllChem as ac
import rdkit
from rdkit.Geometry import Point3D
import re

# This dictionary gives easy access to the rdkit bond types to the 
# ``mol_from_mol2_file`` function.
bond = {'1' : rdkit.Chem.rdchem.BondType.SINGLE,
        'am' : rdkit.Chem.rdchem.BondType.SINGLE,
        '2' : rdkit.Chem.rdchem.BondType.DOUBLE,
        '3' : rdkit.Chem.rdchem.BondType.TRIPLE,
        'ar' : rdkit.Chem.rdchem.BondType.AROMATIC}

def mol_from_mol2_file(mol2_file):
    """
    Creates an rdkit molecule from a ``.mol2`` file.
    
    Parameters
    ----------
    mol2_file : str
        The full path of the ``.mol2`` file from which an rdkit molecule
        should be instantiated.
        
    Returns
    -------
    rdkit.Chem.rdchem.Mol
        An rdkit instance of the molecule held in the `mol2_file`.

    """
    
    # Read the ``.mol2`` file line by line. Checks for the lines
    # holding flags indicating the start of the atomic or bond block.
    # When going through a block use its data in the rdkit molecule or
    # the conformer. Finally add the conformer to the rdkit molecule and
    # return.
    
    mol = chem.Mol()  
    e_mol = chem.EditableMol(mol)
    conf = chem.Conformer()
    atomic_symbol = re.compile('[A-z]')    
    
    
    take_atom = False 
    take_bond = False
    
    with open(mol2_file, 'r') as f:
        for line in f:
            
            # Indicates the following lines hold the atom block.
            if '@<TRIPOS>ATOM' in line:
                take_atom = True
                continue
            # Indicates the following lines hold the bond block, and the
            # atom block as ended.
            if '@<TRIPOS>BOND' in line:
                take_atom  = False
                take_bond = True
                continue
            # Indicates that the bond block is ended and all data has
            # therefore been collected.
            if take_bond and len(line.split()) in {0,1}:
                break
            # If in the atom block, extract atomic data.
            if take_atom:
                _, atom_sym, x, y, z, *_ = line.split()
                atom_sym = atomic_symbol.findall(atom_sym)[0]
                
                atom_id = e_mol.AddAtom(chem.Atom(atom_sym))
                atom_coord = Point3D(float(x), float(y), float(z))                
                
                conf.SetAtomPosition(atom_id, atom_coord)
                
                continue
            
            # If in the bond block, extract bond data.
            if take_bond:
                bond_id, atom1, atom2, bond_order, *_ = line.split()
                e_mol.AddBond(int(atom1)-1, int(atom2)-1, 
                              bond[bond_order])                
                
                continue
    
    # Get the rdkit molecule and give it the conformer.
    mol = e_mol.GetMol()
    mol.AddConformer(conf)
    return mol
    
            
def dedupe(iterable, seen=None):
    if seen is None:
        seen = set()        
    for x in iterable:
        if x not in seen:
            seen.add(x)
            yield x
            
def flatten(iterable, excluded_types={str}):
    for x in iterable:
        if hasattr(x, '__iter__') and type(x) not in excluded_types:          
            yield from flatten(x)
        else:
            yield x            


periodic_table = {1: 'H',
                  2: 'He',
                  3: 'Li',
                  4: 'Be',
                  5: 'B',
                  6: 'C',
                  7: 'N',
                  8: 'O',
                  9: 'F',
                  10: 'Ne',
                  11: 'Na',
                  12: 'Mg',
                  13: 'Al',
                  14: 'Si',
                  15: 'P',
                  16: 'S',
                  17: 'Cl',
                  18: 'Ar',
                  19: 'K',
                  20: 'Ca',
                  21: 'Sc',
                  22: 'Ti',
                  23: 'V',
                  24: 'Cr',
                  25: 'Mn',
                  26: 'Fe',
                  27: 'Co',
                  28: 'Ni',
                  29: 'Cu',
                  30: 'Zn',
                  31: 'Ga',
                  32: 'Ge',
                  33: 'As',
                  34: 'Se',
                  35: 'Br',
                  36: 'Kr',
                  37: 'Rb',
                  38: 'Sr',
                  39: 'Y',
                  40: 'Zr',
                  41: 'Nb',
                  42: 'Mo',
                  43: 'Tc',
                  44: 'Ru',
                  45: 'Rh',
                  46: 'Pd',
                  47: 'Ag',
                  48: 'Cd',
                  49: 'In',
                  50: 'Sn',
                  51: 'Sb',
                  52: 'Te',
                  53: 'I',
                  54: 'Xe',
                  55: 'Cs',
                  56: 'Ba',
                  57: 'La',
                  58: 'Ce',
                  59: 'Pr',
                  60: 'Nd',
                  61: 'Pm',
                  62: 'Sm',
                  63: 'Eu',
                  64: 'Gd',
                  65: 'Tb',
                  66: 'Dy',
                  67: 'Ho',
                  68: 'Er',
                  69: 'Tm',
                  70: 'Yb',
                  71: 'Lu',
                  72: 'Hf',
                  73: 'Ta',
                  74: 'W',
                  75: 'Re',
                  76: 'Os',
                  77: 'Ir',
                  78: 'Pt',
                  79: 'Au',
                  80: 'Hg',
                  81: 'Tl',
                  82: 'Pb',
                  83: 'Bi',
                  84: 'Po',
                  85: 'At',
                  86: 'Rn',
                  87: 'Fr',
                  88: 'Ra',
                  89: 'Ac',
                  90: 'Th',
                  91: 'Pa',
                  92: 'U',
                  93: 'Np',
                  94: 'Pu',
                  95: 'Am',
                  96: 'Cm',
                  97: 'Bk',
                  98: 'Cf',
                  99: 'Es',
                  100: 'Fm',
                  101: 'Md',
                  102: 'No',
                  103: 'Lr',
                  104: 'Rf',
                  105: 'Db',
                  106: 'Sg',
                  107: 'Bh',
                  108: 'Hs',
                  109: 'Mt',
                  110: 'Ds',
                  111: 'Rg',
                  112: 'Cn',
                  113: 'Uut',
                  114: 'Fl',
                  115: 'Uup',
                  116: 'Lv',
                  117: 'Uus',
                  118: 'Uuo'}