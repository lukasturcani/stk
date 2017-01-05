import rdkit.Chem as chem
import rdkit
from rdkit.Geometry import Point3D
import numpy as np
import time
from contextlib import contextmanager
import os
import subprocess as sp 
# Prevents the matplotlib import from printing warnings in Spyder. These
# are printed because Spyder automatically imports matplotlib so the
# call to `matplotlib.use()` results in a warning.
import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

import gzip
import re
from collections import deque
import shutil
from collections import namedtuple
import tarfile

# This dictionary gives easy access to the rdkit bond types.

bond_dict = {'1' : rdkit.Chem.rdchem.BondType.SINGLE,
             'am' : rdkit.Chem.rdchem.BondType.SINGLE,
             '2' : rdkit.Chem.rdchem.BondType.DOUBLE,
             '3' : rdkit.Chem.rdchem.BondType.TRIPLE,
             'ar' : rdkit.Chem.rdchem.BondType.AROMATIC}
        
def archive_output():
    """
    Places the ``output`` folder into ``old_output``.
    
    This function assumes that the ``output`` folder is in the current
    working directory.
    
    Returns
    -------
    None : NoneType
    
    """
    
    if 'output' in os.listdir():
        # Make the ``old_output`` folder if it does not exist already.
        if 'old_output' not in os.listdir():
            os.mkdir('old_output')
        
        # Find out with what number the ``output`` folder should be
        # labelled within ``old_output``.
        num = len(os.listdir('old_output'))
        new_dir = os.path.join('old_output', str(num))
        s = 'Moving old output dir.'
        print('\n'+s + '\n' + '-'*len(s) + '\n\n')
        shutil.copytree('output', new_dir)
    
        # Wait for the copy to complete before removing the old folder.
        mv_complete = False    
        while not mv_complete:
            try:
                shutil.rmtree('output')
                mv_complete = True
            except:
                pass

def tar_output():
    """
    Places all the content in the `outout` into a .gz file.

    This function also deletes all the folders in the `output` folder
    expect the one holding the final generation.

    Returns
    -------
    None : NoneType    
    
    """
    
    s = "Compressing output."
    print("\n\n"+s+"\n"+"-"*len(s)+"\n\n")
    
    with tarfile.open(os.path.join('output','output.gz'), 'w:gz') as tar:
        tar.add('output')
    
    folders = [x for x in os.listdir('output') if 
                os.path.isdir(os.path.join('output', x))]

    max_folder = max((x for x in folders if x != 'initial'), 
                     key=lambda x : int(x))    

    for folder in folders:
        if folder != max_folder:
            shutil.rmtree(os.path.join('output', folder))
        
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

def centroid(*coords):
    """
    Calculates the centroid of a group of coordinates.
    
    Parameters
    ----------
    *coords : numpy.array
        Any number of numpy arrays holding x, y and z positions.
        
    Returns
    -------
    numpy.array
        The centroid of the coordinates `coords`.
    
    """
    
    total = 0
    for coord in coords:
        total = np.add(total, coord)
    return np.divide(total, len(coords))
    

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
    """ 
    
    Returns a rotation matrix of `angle` radians about `axis`.
    
    Parameters
    ----------
    angle : int or float
        The size of the rotation in radians.
        
    axis : numpy.array
        A 3 element aray which represents a vector. The vector is the
        axis about which the rotation is carried out.
        
    Returns
    -------
    numpy.array
        A 3x3 array representing a rotation matrix.
    
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
    matrix : np.matrix
        A n x 3 matrix. Each row holds the x, y and z coordinate of some
        point, respectively.
        
    Returns
    -------
    numpy.array
        A numpy array which holds the x, y and z coordinates of the
        centroid of the coordinates in `matrix`.
    
    """
    

    return np.array(np.sum(matrix, axis=0) / len(matrix))[0]

class ChargedMolError(Exception):
    def __init__(self, mol_file, msg):
        self.mol_file = mol_file
        self.msg = msg
class MolFileError(Exception):
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
    MolFileError
        If the file is not a V3000 .mol file.
    
    """
    
    e_mol = chem.EditableMol(chem.Mol())
    conf = chem.Conformer()
    
    with open(mol_file, 'r') as f:
        take_atom = False
        take_bond = False
        v3000 = False
        
        for line in f:
            if 'V3000' in line:
                v3000 = True
                
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
    if not v3000:
        raise MolFileError(mol_file, 'Not a V3000 .mol file.')
            
    mol = e_mol.GetMol()
    mol.AddConformer(conf)
    return mol

def mol_from_mae_file(mae_path):
    """
    Creates a rdkit molecule from a ``.mae`` file.
    
    Parameters
    ----------
    mol2_file : str
        The full path of the ``.mae`` file from which an rdkit molecule
        should be instantiated.

    Returns
    -------
    rdkit.Chem.rdchem.Mol
        An rdkit instance of the molecule held in `mae_file`.

    """
    
    mol = chem.EditableMol(chem.Mol())
    conf = chem.Conformer()

    with open(mae_path, 'r') as mae:
        content = re.split(r'[{}]', mae.read())
    
    prev_block = deque([''], maxlen=1)
    for block in content:
        if 'm_atom[' in prev_block[0]:
            atom_block = block
        if 'm_bond[' in prev_block[0]:
            bond_block = block
        prev_block.append(block)
    


    labels, data_block, *_ = atom_block.split(':::')
    labels = [l for l in labels.split('\n') if 
               not l.isspace() and l != '']
    
    data_block = [a.split() for a in data_block.split('\n') if 
                   not a.isspace() and a != '']    
    
    for line in data_block:
        line = [word for word in line if word != '"']
        if len(labels) != len(line):
            raise RuntimeError(('Number of labels does'
                      ' not match number of columns in .mae file.'))
            
        for label, data in zip(labels, line):
            if 'x_coord' in label:
                x = float(data)
            if 'y_coord' in label:
                y = float(data)
            if 'z_coord' in label:
                z = float(data)
            if 'atomic_number' in label:
                atom_num = int(data)
        
        atom_sym = periodic_table[atom_num]        
        atom_coord = Point3D(x,y,z)
        atom_id = mol.AddAtom(chem.Atom(atom_sym))              
        conf.SetAtomPosition(atom_id, atom_coord)        


    
    labels, data_block, *_ = bond_block.split(':::')
    labels = [l for l in labels.split('\n') if not l.isspace() and l != '']
    data_block = [a.split() for a in data_block.split('\n') if not a.isspace() and a != '']
    
    for line in data_block:
        if len(labels) != len(line):
            raise RuntimeError(('Number of labels does'
                      ' not match number of columns in .mae file.'))
            
        for label, data in zip(labels, line):
            if 'from' in label:
                atom1 = int(data) - 1
            if 'to' in label:
                atom2 = int(data) - 1
            if 'order' in label:
                bond_order = str(int(data))
        mol.AddBond(atom1, atom2, bond_dict[bond_order])  
        
    mol = mol.GetMol()
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
    x_vals = list(range(1, len(counter.items())+1))
    y_vals = []
    labels = []
    
    for ind, value in sorted(counter.items(), reverse=True):
        y_vals.append(value)
        labels.append(ind.fitness)
    
    plt.bar(x_vals, y_vals, color='blue')
    plt.xticks([x+0.5 for x in x_vals], labels, rotation='vertical')
    plt.tight_layout()
    fig.savefig(plot_name, dpi=fig.dpi)
    plt.close('all')

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
    print('\nTime taken was {0} : {1} : {2}.\n\n'.format(
                                                    int(h), int(m), s))


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
       
class MAEExtractor:
    """
    Extracts the lowest energy conformer from a .maegz file.

    Attributes
    ----------
    macro_mol : MacroMolecule 
            
    """
    
    def __init__(self, macro_mol):
        self.maegz_path = macro_mol.prist_mol_file.replace('.mol', 
                                                           "-out.maegz")
        self.macro_mol = macro_mol
        self.maegz_to_mae()
        self.extract_conformer()

    def extract_conformer(self):
        print("Extracting conformer - {}.".format(
                                        self.macro_mol.prist_mol_file))        
                
        
        num = self.lowest_energy_conformer()        
        
        content = self.content.split("f_m_ct")
        new_mae = "f_m_ct".join([content[0], content[num]])
        new_name = self.mae_path.replace('.mae', '_EXTRACTED_{}.mae'.format(num))
        with open(new_name, 'w') as mae_file:
            mae_file.write(new_mae)
        self.path = new_name
            
    def extract_energy(self, block):
        block = block.split(":::")
        for name, value in zip(block[0].split('\n'), block[1].split('\n')):
            if 'r_mmod_Potential_Energy' in name:
                return float(value)            
        
    def lowest_energy_conformer(self):
        with open(self.mae_path, 'r') as mae_file:
            self.content = mae_file.read()
            self.content_split = re.split(r"[{}]", self.content)

        energies = []

        prev_block = deque([""], maxlen=1)
        index = 1        
        for block in self.content_split:
            if ("f_m_ct" in prev_block[0] and
                                "r_mmod_Potential_Energy" in block):              
                energy = self.extract_energy(block)
                energies.append((energy, index))
                index += 1
                
            prev_block.append(block)
        return min(energies)[1]
  
    def maegz_to_mae(self):
        self.mae_path = self.maegz_path.replace('.maegz', '.mae')            
        with gzip.open(self.maegz_path, 'r') as maegz_file:
            with open(self.mae_path, 'wb') as mae_file:
                mae_file.write(maegz_file.read())       
                
class FGInfo:
    """
    Contains key information for functional group substitutions.
    
    The point of this class is to register which atom is substituted
    for which, when an atom in a functional group is substituted with a 
    heavy metal atom. If MMEA is to incorporate a new functional group, 
    a new ``FGInfo`` instance should be added to the 
    `functional_group_list` class attribute of ``FGInfo``. 
    
    Adding a new ``FGInfo`` instace to `functional_group_list` will 
    allow the `Topology.join_mols` method to connect this functional 
    group to (all) others during assembly. Nothing except adding this
    instance should need to be done in order to incorporate new 
    functional groups.
    
    If this new functional group is to connect to another functional 
    group with a double bond during assembly, the symbols of the heavy 
    atoms of both functional groups should be added to the 
    `double_bond_combs` class attribute. The order in which the heavy 
    symbols are placed in the tuple does not matter. Again, this is all
    that needs to be done for MMEA to create double bonds between
    certain functional groups.  
    
    Class attributes
    ----------------
    functional_groups_list : list of FGInfo instances
        This list holds all ``FGInfo`` instances used by MMEA. If a new
        functional group is to be used by MMEA, a new ``FGInfo`` 
        instance must be added to this list.
        
    double_bond_combs : list of tuples of strings
        When assembly is carried out, if the heavy atoms being joined
        form a tuple in this list, they will be joined with a double
        rather than single bond. If a single bond is desired there is no
        need to change this variable.
        
    heavy_symbols : set of str
        A set of all the heavy symbols used by ``FGInfo`` instances in 
        MMEA. This set updates itself automatically. There is no need to
        modify it when changes are made to any part of MMEA.
        
    heavy_atomic_nums : set of ints
        A set of all atomic numbers of heavy atoms used by ``FGInfo``
        instances in MMEA. This set updates itself automatically. There
        is no need to modify it when chagnes are made to any part of
        MMEA.

    Attributes
    ----------
    name : str
        The name of the functional group.
    
    smarts_start : str
        A ``SMARTS`` string describing the functional group before 
        substitution by a heavy atom.
        
    del_tags : list of DelAtom instances
        Every member of this list represents an atom on the functional
        group which should be deleted during assembly. One atom in each
        functional group is removed for each list member.
    
    target_atomic_num : int
        The atomic number of the atom being substituted by a heavy atom.
    
    heavy_atomic_num : int
        The atomic number of the heavy atom which replaces the target 
        atom.
    
    target_symbol : str
        The atomic symbol of the atom, being substituted by a heavy 
        atom.       
    
    heavy_symbol : str
        The atomic symbol of the heavy atom which replaces the target 
        atom.
    
    """
    
    __slots__ = ['name', 'smarts_start', 'del_tags', 
                 'target_atomic_num', 'heavy_atomic_num', 
                 'target_symbol', 'heavy_symbol'] 
    
    def __init__(self, name, smarts_start, del_tags, target_atomic_num, 
                 heavy_atomic_num, target_symbol, heavy_symbol):
         self.name = name
         self.smarts_start = smarts_start
         self.del_tags = del_tags
         self.target_atomic_num = target_atomic_num
         self.heavy_atomic_num = heavy_atomic_num
         self.target_symbol = target_symbol
         self.heavy_symbol = heavy_symbol

# An atom is deleted based on what type of bond connects it to the
# substituted functional group atom. The element of the atom is ofcourse
# a factor as well. When both of these are satisfied the atom is
# removed. The ``DelAtom`` class conveniently stores this information.
# Bond type is an rdkit bond type (see the bond dictionary above for
# the two possible values it may take) and atomic num in an integer.
DelAtom = namedtuple('DelAtom', ['bond_type', 'atomic_num'])

FGInfo.functional_group_list = [
                        
    FGInfo("aldehyde", "C(=O)[H]", [ DelAtom(bond_dict['2'], 8) ], 
                                                       6, 39, "C", "Y"), 
    
    FGInfo("carboxylic_acid", "C(=O)O[H]", 
           [ DelAtom(bond_dict['1'], 8) ], 6, 40, "C", "Zr"),
    
    FGInfo("amide", "C(=O)N([H])[H]", [ DelAtom(bond_dict['1'], 7) ], 
                                                      6, 41, "C", "Nb"),
    
    FGInfo("thioacid", "C(=O)S[H]", [ DelAtom(bond_dict['1'], 16) ], 
                                                      6, 42, "C", "Mo"),
    
    FGInfo("alcohol", "O[H]", [], 8, 43, "O", "Tc"),
    FGInfo("thiol", "[S][H]", [], 16, 44, "S", "Ru"),
    FGInfo("amine", "[N]([H])[H]", [], 7, 45, "N", "Rh"),       
    FGInfo("nitroso", "N=O", [], 7, 46, "N", "Pd"),
    FGInfo("boronic_acid", "[B](O[H])O[H]", [], 5, 47, "B", "Ag")
                             
                             ]

FGInfo.double_bond_combs = [("Rh","Y"), ("Nb","Y"), ("Mb","Rh")]

FGInfo.heavy_symbols = {x.heavy_symbol for x 
                                        in FGInfo.functional_group_list}
                        
FGInfo.heavy_atomic_nums = {x.heavy_atomic_num for x 
                                        in FGInfo.functional_group_list}
def kill_macromodel():
    """
    Kills any applications left open as a result running MacroModel.    
    
    Applications that are typically left open are 
    ``jserver-watcher.exe`` and ``jservergo.exe``.    
    
    Returns
    -------
    None : NoneType    
    
    """
    
    if os.name == 'nt':
        # In Windows, use the ``Taskkill`` command to force a close on
        # the applications.           
        sp.run(["Taskkill", "/IM", "jserver-watcher.exe", "/F"], 
               stdout=sp.PIPE, stderr=sp.PIPE)
        sp.run(["Taskkill", "/IM", "jservergo.exe", "/F"],
               stdout=sp.PIPE, stderr=sp.PIPE)
    if os.name == 'posix':
        sp.run(["pkill", "jservergo"],
               stdout=sp.PIPE, stderr=sp.PIPE)
        sp.run(["pkill", "jserver-watcher"],
               stdout=sp.PIPE, stderr=sp.PIPE)    
                
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
