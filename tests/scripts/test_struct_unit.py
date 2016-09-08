import os
import rdkit
import numpy as np
import copy

from ...classes import StructUnit, Linker, FGInfo
from ...convenience_functions import flatten, normalize_vector

def get_mol_file():
    # The following lines first create a directory tree starting from
    # the current working directory. A generator expression then
    # iterates through the directory, where the ``if`` condition ensures
    # that the desired ``.mol`` file is found. It is the one in the 
    # ``data`` directory. Finally the full path of the ``.mol`` file is 
    # generated using the ``os.path.join`` function. This approach means
    # that the test should work on any machine as it does not depend on
    # absolute paths to find the ``.mol`` file. It also means that the 
    # test does not need to be run from a specific directory.    
    mol_file = os.walk(os.getcwd())    
    for x in mol_file:
        if 'data' in x[0]:
            for y in x[2]:
                if '.mol' in y and 'HEAVY' not in y:
                    yield os.path.join(x[0], y)

# Create a StructUnit instance which can be used by multiple tests to
# not waste time on multiple initializations.
struct_file = next(x for x in get_mol_file() if 'amine3f_14.mol' in x)
struct = StructUnit(struct_file)

def test_caching():
    bb_file = next(x for x in get_mol_file() 
                                    if 'amine3f_14.mol' in x)
    lk_file = next(x for x in get_mol_file() 
                                    if 'aldehyde2f_3.mol' in x) 

    bb = StructUnit(bb_file)
    bb2 = StructUnit(bb_file)
    
    lk = StructUnit(lk_file)
    lk2 = StructUnit(lk_file)
    
    assert bb is bb2
    assert bb is not lk
    assert lk is lk2
    assert lk2 is not bb2

def test_init():
    """
    Ensures that StructUnit instances are initiated correctly.
    
    This function only checks that the attributes exist and that
    their values of of the correct type. Checking that correct values
    are initialized is done in other tests.

    """     
    mol_file = next(x for x in get_mol_file() 
                                        if 'aldehyde2f_3.mol' in x)
    struct_unit = StructUnit(mol_file)
     
    # Check that heavy attributes were created by the initializer.
    assert hasattr(struct_unit, 'heavy_mol')
    assert hasattr(struct_unit, 'heavy_mol_file')
    assert hasattr(struct_unit, 'heavy_smiles')
    
    assert isinstance(struct_unit.func_grp, FGInfo)
    assert isinstance(struct_unit.prist_mol, rdkit.Chem.rdchem.Mol)
    assert isinstance(struct_unit.heavy_mol, rdkit.Chem.rdchem.Mol)
    assert isinstance(struct_unit.prist_smiles, str)
    assert isinstance(struct_unit.heavy_smiles, str)
    assert isinstance(struct_unit.prist_mol_file, str)
    assert isinstance(struct_unit.heavy_mol_file, str)
    
def test_find_functional_group_atoms():
    """
    Make sure correct atoms are found in the functional groups.
    
    This function uses the ``aldehyde2f_3.mol`` file.

    """    
    # These are the expected atom ids for the test molecule.
    expected = ((1, 0, 12), (10, 11, 19))    
    
    # Initializing the test molecule.    
    mol_file = next(x for x in get_mol_file() 
                                        if 'aldehyde2f_3.mol' in x)
    struct_unit = StructUnit(mol_file)
        
    func_grp_atoms = struct_unit.find_functional_group_atoms()

    assert func_grp_atoms == expected
    
def test_shift_heavy_mol():
    """
    Ensure that shifting the position of a ``StructUnit`` works.    
    
    This function uses the ``aldehyde2f_3.mol`` file.    
    
    """
    # Initializing the test molecule.    
    mol_file = next(x for x in get_mol_file() 
                                        if 'aldehyde2f_3.mol' in x)
    struct_unit = StructUnit(mol_file)

    # Shifting the same molecule twice should return two 
    # ``rdkit.Chem.rdchem.Mol`` instances with conformers describing the
    # same atomic positions. Furthermore, the original conformer should
    # in the ``StructUnit`` instance should be unchanged.
    og_conformer = struct_unit.heavy_mol.GetConformer()

    shift_size = 10    
    a = struct_unit.shift_heavy_mol(shift_size, shift_size, shift_size)
    a_conformer = a.GetConformer()
    
    b = struct_unit.shift_heavy_mol(shift_size, shift_size, shift_size)
    b_conformer = b.GetConformer()
    
    # Check that the same atom coords are present in `a` and `b`. Also
    # Check that these are different to the coords in the original.
    for atom in a.GetAtoms():
        atom_id = atom.GetIdx()
      
        og_coords = og_conformer.GetAtomPosition(atom_id)
        atom1_coords = a_conformer.GetAtomPosition(atom_id)
        atom2_coords = b_conformer.GetAtomPosition(atom_id)
        
        assert atom1_coords.x == atom2_coords.x
        assert atom1_coords.y == atom2_coords.y
        assert atom1_coords.z == atom2_coords.z
        
        # Checking that the coords are differnt to original. By the 
        # correct amount.        
        assert atom1_coords.x == og_coords.x + shift_size
        assert atom1_coords.y == og_coords.y + shift_size
        assert atom1_coords.z == og_coords.z + shift_size
        
def test_heavy_all_atom_coords():
    """
    Make sure the correct output is provided.

    """
    mol_file = next(x for x in get_mol_file() 
                                        if 'amine3f_14.mol' in x)   
    mol = StructUnit(mol_file)
    expected_output_type = type(np.array([1,2,3]))
    
    for atom_id, coord in mol.heavy_all_atom_coords():
        x,y,z = coord
        
        assert isinstance(atom_id, int)
        assert isinstance(coord, expected_output_type)        
        assert isinstance(x, float)
        assert isinstance(y, float)
        assert isinstance(z, float)
    
    assert len(list(mol.heavy_all_atom_coords())) == 32
    
def test_amine_substitution():
    """
    Ensure that the amine functional group is correctly replaced.    
    
    """
    exp_smiles = ("[H][C]([H])([H])[C]1=[C]([Rh])[C](=[O])[C]2=[C]"
                  "([C]1=[O])[N]1[C](=[C]2[C]([H])([H])[O][C](=[O])"
                  "[Rh])[C]([H])([H])[C]([H])([Rh])[C]1([H])[H]")
    mol_file = next(x for x in get_mol_file() 
                                        if 'amine3f_14.mol' in x)   
    mol = StructUnit(mol_file)
    assert mol.heavy_smiles == exp_smiles

def test_aldehyde_substitution():
    """
    Ensure that the aldehyde functional group is correctly replaced.    
    
    """
    exp_smiles = ("[H][N]1/[C](=[N]/[Y])[C]([H])([H])[N]([H])[C]([H])"
                    "([H])/[C]1=[N]\[Y]")
    mol_file = next(x for x in get_mol_file() 
                                        if 'aldehyde2f_3.mol' in x)   
    mol = StructUnit(mol_file)
    print(mol.heavy_smiles)
    assert mol.heavy_smiles == exp_smiles
    
def test_make_atoms_heavy_in_heavy():
    """
    This test might need more assert statements.
    
    """
    

    lk_file = next(x for x in get_mol_file() 
                                    if 'aldehyde2f_3.mol' in x) 

    lk = StructUnit(lk_file)
    
    # Test that the position of the substituted atoms remains the same.
        
    i= 0
    for atom_id in flatten(lk.find_functional_group_atoms()):
        atom = lk.prist_mol.GetAtomWithIdx(atom_id)
        if atom.GetAtomicNum() == lk.func_grp.target_atomic_num:
            prist_coord = lk.prist_get_atom_coords(atom_id)
            heavy_coord = lk.heavy_get_atom_coords(lk.heavy_ids[i])
            i += 1
            assert np.array_equal(prist_coord, heavy_coord)        
    
def test_set_heavy_mol_position():

        
    # Place centroid at some position.
    mol = struct.set_heavy_mol_position([3.14, 6.14, -12.14])
    
    # Check that the centroid of the heavy molecule is at that position.
    assert np.allclose(struct.heavy_mol_centroid(), 
                       np.array([3.14, 6.14, -12.14]),
                       atol = 1e-8)
    
    # Ensure the returned rdkit instance is the one in `heavy_mol`.                   
    assert mol is struct.heavy_mol

def test_heavy_mol_position_matrix():
    
    # Go through each atom id. For each atom id get the column in the 
    # position matrix with that id as its index. Make sure that the data
    # is the same. 
    pos_mat1 = struct.heavy_mol_position_matrix()
    conf = struct.heavy_mol.GetConformer()
       
    for atom in struct.heavy_mol.GetAtoms():
        atom_id = atom.GetIdx()
        cx, cy, cz = conf.GetAtomPosition(atom_id)
        
        conf_coord = np.array([cx, cy, cz])   
        mat_coord = pos_mat1.T[atom_id]

        assert np.allclose(conf_coord, mat_coord, atol = 1e-8)
    
    # Move the molecule, ensure that the position matrix is adjusted
    # appropriately.
    curr_x, curr_y, curr_z = struct.heavy_mol_centroid()
    
    struct.set_heavy_mol_position([curr_x+1, curr_y-1, curr_z+2])
    pos_mat2 = struct.heavy_mol_position_matrix() 
    conf = struct.heavy_mol.GetConformer()

    for atom in struct.heavy_mol.GetAtoms():
        atom_id = atom.GetIdx()
        cx, cy, cz = conf.GetAtomPosition(atom_id)
        
        conf_coord = np.array([cx, cy, cz])   
        mat_coord = pos_mat1.T[atom_id]
        mat2_coord = pos_mat2.T[atom_id]
        
        assert not np.allclose(conf_coord, mat_coord, atol = 1e-8)    
        assert np.allclose(conf_coord, mat2_coord, atol = 1e-8) 
    
def test_set_heavy_mol_from_position_matrix_AND_OTHERS():
    """
    Also tests heavy_mol_centroid and heavy_atom_centroid.
    
    """
    
    # Make a position matrix where each coordinate is set to [1,2,3] for
    # every atom. Use a copy of struct in this test because it will
    # severly mess up the structure.
    
    struct2 = copy.deepcopy(struct)
    
    pos_mat = []
    for x in range(struct2.heavy_mol.GetNumAtoms()):
        pos_mat.append([1,2,3])
    
    pos_mat = np.matrix(pos_mat).T

    struct2.set_heavy_mol_from_position_matrix(pos_mat)
    
    for _, coord in struct2.heavy_all_atom_coords():
        assert np.array_equal(coord, [1,2,3])
    
    # Centroids should also be in [1,2,3].
    assert np.array_equal(struct2.heavy_atom_centroid(), [1,2,3])
    assert np.array_equal(struct2.heavy_mol_centroid(), [1,2,3])

def test_set_heavy_mol_orientation():
    
    struct.set_heavy_mol_orientation(
    next(struct.heavy_direction_vectors()), [1,2,3])
    
    assert np.allclose(next(struct.heavy_direction_vectors()), 
                       normalize_vector([1,2,3]),
                       atol=1e-8)   
    
def test_heavy_atom_position_matrix():
    
    conf = struct.heavy_mol.GetConformer()
    pos_mat = struct.heavy_atom_position_matrix()    
    
    for atom_id in struct.heavy_ids:
        cx, cy, cz = conf.GetAtomPosition(atom_id)
        coord = np.array([cx, cy, cz])        
        
        column_i = struct.heavy_ids.index(atom_id)
        column = pos_mat.T[column_i]
        
        assert np.allclose(coord, column, atol=1e-8)
        

def test_heavy_direction_vectors():
    
    for vector in struct.heavy_direction_vectors():
        assert isinstance(vector, np.ndarray)
    
    assert len(list(struct.heavy_direction_vectors())) == 3
    
    
def test_set_heavy_atom_centroid():
    struct.set_heavy_atom_centroid([1,2,3])
    assert np.allclose(struct.heavy_atom_centroid(), [1,2,3], atol=1e-8)
    
    struct.set_heavy_atom_centroid([2,2,2])
    assert not np.allclose(struct.heavy_atom_centroid(), 
                           [1,2,3], atol=1e-8)     

    assert np.allclose(struct.heavy_atom_centroid(), [2,2,2], atol=1e-8)    
    
def test_prist_get_atom_coords():
    
    conf = struct.prist_mol.GetConformer()    
    
    for atom in struct.prist_mol.GetAtoms():
        atom_id = atom.GetIdx()
        coords = struct.prist_get_atom_coords(atom_id)
        conf_coords = conf.GetAtomPosition(atom_id)
        assert np.allclose(coords, conf_coords, atol=1e-8)
    
def test_heavy_get_atom_coords():
    
    conf = struct.heavy_mol.GetConformer()    
    
    for atom in struct.heavy_mol.GetAtoms():
        atom_id = atom.GetIdx()
        coords = struct.heavy_get_atom_coords(atom_id)
        conf_coords = conf.GetAtomPosition(atom_id)
        assert np.allclose(coords, conf_coords, atol=1e-8)
    
    
    
    