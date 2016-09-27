import rdkit.Chem as chem
import rdkit
from rdkit.Geometry import Point3D
import numpy as np
import time
from contextlib import contextmanager
import matplotlib.pyplot as plt
import os 

# This dictionary gives easy access to the rdkit bond types.
bond_dict = {'1' : rdkit.Chem.rdchem.BondType.SINGLE,
             'am' : rdkit.Chem.rdchem.BondType.SINGLE,
             '2' : rdkit.Chem.rdchem.BondType.DOUBLE,
             '3' : rdkit.Chem.rdchem.BondType.TRIPLE,
             'ar' : rdkit.Chem.rdchem.BondType.AROMATIC}

def normalize_vector(vector):
    """
    Normalizes the given vector.
    
    A new vector is returned, the original vector is not modified.    
    
    Parameters
    ----------
    vector : np.array
        The vector to be normalized.        
        
    Returns
    -------
    np.array
        The normalized vector.
    
    """

    v = np.divide(vector, np.linalg.norm(vector))
    return np.round(v, decimals=4)

def vector_theta(vector1, vector2):
    """
    Returns the angle between two vectors in radians.

    Parameters
    ----------
    vector1 : numpy.array
        The first vector.
        
    vector2 : numpy.array    
        The second vector.
        
    Returns
    -------
    float
        The angle between `vector1` and `vector2` in radians.
    
    """
    
    numerator = np.dot(vector1, vector2)
    denominator = (np.linalg.norm(vector1) * 
                    np.linalg.norm(vector2))
    
    return np.arccos(numerator/denominator)    

def rotation_matrix(vector1, vector2):
    """
    Returns a rotation matrix which transforms `vector1` to `vector2`.

    Multiplying `vector1` by the rotation matrix returned by this 
    function yields `vector2`. 

    Parameters
    ----------
    vector1 : numpy.array
        The vector which needs to be transformed to `vector2`.

    vector2 : numpy.array
        The vector onto which `vector1` needs to be transformed.
    
    Returns
    -------
    numpy.ndarray
        A rotation matrix which transforms `vector1` to `vector2`.
        
    References
    ----------
    http://tinyurl.com/kybj9ox
    http://tinyurl.com/gn6e8mz
    
    """

    # Make sure both inputs are unit vectors.
    vector1 = normalize_vector(vector1)
    vector2 = normalize_vector(vector2)
    
    # Hande the case where the input and output vectors are equal.
    if np.array_equal(vector1, vector2):
        return np.identity(3)
    
    # Handle the case where the rotation is 180 degrees.
    if np.array_equal(vector1, np.multiply(vector2, -1)):
        return np.multiply(np.identity(3), -1)
        
    v = np.cross(vector1, vector2)
    
    vx = np.array([[0, -v[2], v[1]],
                   [v[2], 0, -v[0]], 
                   [-v[1], v[0], 0]])
    
    s = np.linalg.norm(v)
    c = np.dot(vector1, vector2)
    I = np.identity(3)

    mult_factor = (1-c)/np.square(s)    
    
    return I + vx + np.multiply(np.dot(vx,vx), mult_factor)

def kabsch(coords1, coords2):
    """
    Return a rotation matrix to minimize dstance between 2 coord sets.
    
    This is essentially an implementation of the Kabsch algorithm. Given
    two sets of coordinates, `coords1` and `coords2`, this function
    returns a rotation matrix. When the rotation matrix is applied to
    `coords1` the resulting coordinates have their rms distance to
    `coords2` minimized.
    
    Parameters
    ----------
    coords1 : numpy.array
        This array represents a matrix hodling coordinates which need
        to be rotated to minimize their rms distance to coordinates in
        `coords2`. The matrix is n x 3. Each row of the matrix holds the
        x, y and z coordinates of one point, respectively. Here ``n`` is
        the number of points.
    
    coords2 : numpy.array
        This array represents a matrix which holds the coordinates of
        points the distance to which should be minimized. The matrix is 
        n x 3. Each row of the matrix holds the x, y and z coordinates 
        of one point, respectively. Here ``n`` is the number of points.
        
    Returns
    -------
    numpy.array
        A rotation matrix. This will be a 3 x 3 matrix.
    
    References
    ----------
    http://nghiaho.com/?page_id=671
    https://en.wikipedia.org/wiki/Kabsch_algorithm
    
    """
    
    h = np.dot(coords1, coords2.T)
    u,s,v = np.linalg.svd(h)

    if int(np.linalg.det(v)) < 0:
        v[:,2] = -v[:,2]

    
    return np.dot(v, u)

def rotation_matrix_arbitrary_axis(angle, axis):
    """ Applies the rotation to a molecule by using the quaternion method. 
    The coordinates of each atom are rotated by a specific Angle is while 
    keeping the axis fixed. It return the rotated coordinates for the molecule.
    """
    # Calculation of the rotation matrix
    axis = normalize_vector(axis)    
    
    a = np.cos(angle/2)

    b,c,d = np.multiply(axis, np.sin(angle/2))
    
    e11 = np.square(a) + np.square(b) - np.square(c) - np.square(d)    
    e12 = 2*(np.multiply(b,c) - np.multiply(a,d))
    e13 = 2*(np.multiply(b,d) + np.multiply(a,c))
    
    e21 = 2*(np.multiply(b,c) + np.multiply(a,d))
    e22 = np.square(a) + np.square(c) - np.square(b) - np.square(d)
    e23 = 2*(np.multiply(c,d) - np.multiply(a,b))
    
    e31 = 2*(np.multiply(b,d) - np.multiply(a,c))
    e32 =  2*(np.multiply(c,d) + np.multiply(a,b))
    e33 = np.square(a) + np.square(d) - np.square(b) - np.square(c)
    
    return np.array([[e11, e12, e13],
                     [e21, e22, e23],
                     [e31, e32, e33]]) 

    
def matrix_centroid(matrix):
    """
    Returns the centroid of the coordinates held in `matrix`.
    
    Parameters
    ----------
    matrix : np.array
        A n x 3 matrix. Each row holds the x, y and z coordinate of some
        point, respectively.
        
    Returns
    -------
    numpy.array
        A numpy array which holds the x, y and z coordinates of the
        centroid of the coordinates in `matrix`.
    
    """
    
    sum_ = sum( x[0] for x in matrix)
    return sum_ / len(matrix)

class ChargedMolError(Exception):
    def __init__(self, mol_file, msg):
        self.mol_file = mol_file
        self.msg = msg
    

def mol_from_mol_file(mol_file):
    """
    Creates a rdkit molecule from a ``.mol`` (V3000) file.
    
    Parameters
    ----------
    mol_file : str
        The full of the .mol file from which an rdkit molecule should
        be instantiated.
        
    Returns
    -------
    rdkit.Chem.rdchem.Mol
        An rdkit instance of the molecule held in `mol2_file`.
        
    Raises
    ------
    ChargedMolError
        If an atom row has more than 8 coloumns it is usually because
        there is a 9th coloumn indicating atomic charge. Such molecules
        are not currently supported, so an error is raised.
    
    """
    
    mol = chem.Mol()  
    e_mol = chem.EditableMol(mol)
    conf = chem.Conformer()
    
    with open(mol_file, 'r') as f:
        take_atom = False
        take_bond = False
        
        for line in f:
            if 'M  V30 BEGIN ATOM' in line:
                take_atom = True
                continue
            
            if 'M  V30 END ATOM' in line:
                take_atom = False
                continue

            if 'M  V30 BEGIN BOND' in line:
                take_bond = True
                continue
            
            if 'M  V30 END BOND' in line:
                take_bond = False
                continue
            
            if take_atom:
                words = line.split()
                if len(words) > 8:
                    raise ChargedMolError(mol_file, ('Atom row has more'
                    ' than 8 coloumns. Likely due to a charged atom.'))
                _, _, _, atom_sym, *coords, _ = words
                coords = [float(x) for x in coords]
                atom_coord = Point3D(*coords)   
                atom_id = e_mol.AddAtom(chem.Atom(atom_sym))              
                conf.SetAtomPosition(atom_id, atom_coord)
                continue
            
            if take_bond:
                *_, bond_id,  bond_order, atom1, atom2 = line.split()
                e_mol.AddBond(int(atom1)-1, int(atom2)-1, 
                              bond_dict[bond_order])                   
                continue

    mol = e_mol.GetMol()
    mol.AddConformer(conf)
    return mol

def mol_from_mol2_file(mol2_file):
    """
    Creates a rdkit molecule from a ``.mol2`` file.
    
    Parameters
    ----------
    mol2_file : str
        The full path of the ``.mol2`` file from which an rdkit molecule
        should be instantiated.

    Returns
    -------
    rdkit.Chem.rdchem.Mol
        An rdkit instance of the molecule held in `mol2_file`.

    """
    
    # Read the ``.mol2`` file line by line. Checks for the lines
    # holding flags indicating the start of the atomic or bond block.
    # When going through a block use its data in the rdkit molecule or
    # the conformer. Finally add the conformer to the rdkit molecule and
    # return.
    
    # Create an empty molecule and make it editable.
    mol = chem.Mol()  
    e_mol = chem.EditableMol(mol)
    # Create a new conformer.
    conf = chem.Conformer()
      
    take_atom = False 
    take_bond = False
    
    with open(mol2_file, 'r') as f:
        for line in f:
            
            # Indicates the following lines hold the atom block.
            if '@<TRIPOS>ATOM' in line:
                take_atom = True
                continue
            # Indicates the following lines hold the bond block, and the
            # atom block has ended.
            if '@<TRIPOS>BOND' in line:
                take_atom  = False
                take_bond = True
                continue
            # Indicates that the bond block is ended and all data has
            # therefore been collected. Stop going through the file as a
            # result.
            if take_bond and len(line.split()) in {0,1}:
                break
            # If in the atom block, extract atomic data.
            if take_atom:
                _, atom_sym, x, y, z, *_ = line.split()
                atom_sym = ''.join(x for x in atom_sym if x.isalpha())
                atom_id = e_mol.AddAtom(chem.Atom(atom_sym))
                atom_coord = Point3D(float(x), float(y), float(z))                
                
                conf.SetAtomPosition(atom_id, atom_coord)
                
                continue
            
            # If in the bond block, extract bond data.
            if take_bond:
                bond_id, atom1, atom2, bond_order, *_ = line.split()
                e_mol.AddBond(int(atom1)-1, int(atom2)-1, 
                              bond_dict[bond_order])                
                
                continue
    
    # Get the rdkit molecule and give it the conformer.
    mol = e_mol.GetMol()
    mol.AddConformer(conf)
    return mol

def del_charged_mols(database):
    """
    Removes charged .mol files from a folder of .mol files.
    
    Parameters
    ----------
    database : str
        The full path of the folder containing the .mol files.
    
    Returns
    -------
    None : NoneType    
    
    """
    
    for mol_file in os.listdir(database):
        full_path = os.path.join(database, mol_file)
        with open(full_path, 'r') as f:
            content = f.read()
        if 'CHG=' in content:
            print('Deleting', mol_file)
            os.remove(full_path)
    
def plot_counter(counter, plot_name):
    """
    Saves a .png file holding a plot of `counter`.
    
    The counter should hold the number of times a certain population
    member was selected.
    
    Parameters
    ----------
    counter : Counter
        A counter of which members of a population were picked.
    
    plot_name : str
        The full path of the .png where the plot is to be saved.
        
    Returns
    -------
    None : NoneType
    
    """
    fig = plt.figure()
    x_vals = [x.fitness for x in counter.keys()]
    y_vals = list(counter.values())
    plt.bar(x_vals, y_vals)
    fig.savefig(plot_name, dpi=fig.dpi)

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

@contextmanager
def time_it():
    start = time.time()  
    yield
    time_taken = time.time() - start
    m,s = divmod(time_taken, 60)
    h,m = divmod(m, 60)
    print('\nTime taken was {0} : {1} : {2}'.format(int(h), int(m), s))


class LazyAttr:
    """
    A descriptor for creating lazy attributes.
    
    """
    
    def __init__(self, func):
        self.func = func
        
    def __get__(self, obj, cls):
        if obj is None:
            return self
        val = self.func(obj)
        setattr(obj, self.func.__name__, val)
        return val
        
        
# A dictionary which matches atomic number to elemental symbols.
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
                  
atom_vdw_radii = {
                  'Al': 2,    'Sb': 2,    'Ar': 1.88, 'As': 1.85, 'Ba': 2,    'Be': 2,    'Bi': 2, 
                  'B':  2,    'Br': 1.85, 'Cd': 1.58, 'Cs': 2,    'Ca': 2,    'C':  1.7,  'Ce': 2, 
                  'Cl': 1.75, 'Cr': 2,    'Co': 2,    'Cu': 1.4,  'Dy': 2,    'Er': 2,    'Eu': 2,
                  'F':  1.47, 'Gd': 2,    'Ga': 1.87, 'Ge': 2,    'Au': 1.66, 'Hf': 2,    'He': 1.4,
                  'Ho': 2,    'H':  1.09, 'In': 1.93, 'I':  1.98, 'Ir': 2,    'Fe': 2,    'Kr': 2.02,
                  'La': 2,    'Pb': 2.02, 'Li': 1.82, 'Lu': 2,    'Mg': 1.73, 'Mn': 2,    'Hg': 1.55,
                  'Mo': 2,    'Nd': 2,    'Ne': 1.54, 'Ni': 1.63, 'Nb': 2,    'N':  1.55, 'Os': 2,
                  'O':  1.52, 'Pd': 1.63, 'P':  1.8,  'Pt': 1.72, 'K':  2.75, 'Pr': 2,    'Pa': 2,
                  'Re': 2,    'Rh': 2,    'Rb': 2,    'Ru': 2,    'Sm': 2,    'Sc': 2,    'Se': 1.9,
                  'Si': 2.1,  'Ag': 1.72, 'Na': 2.27, 'Sr': 2,    'S':  1.8,  'Ta': 2,    'Te': 2.06,
                  'Tb': 2,    'Tl': 1.96, 'Th': 2,    'Tm': 2,    'Sn': 2.17, 'Ti': 2,    'W':  2,
                  'U':  1.86, 'V':  2,    'Xe': 2.16, 'Yb': 2,    'Y':  2,    'Zn': 1.29, 'Zr': 2,
                  'X':  1.0,  'D':  1.0
                 }                  