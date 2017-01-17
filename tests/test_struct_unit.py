from os.path import join
import rdkit
import numpy as np
import pickle
import pytest
import itertools as it
from scipy.spatial.distance import euclidean

from ...classes import StructUnit, FGInfo
from ...convenience_functions import (flatten, normalize_vector,
                                     vector_theta)

data_dir = join('data', 'molecule', 'molecule.mol')
mol = StructUnit(data_dir)

def test_bonder_atom_position_matrix():
    """
    Tests the output of `bonder_atom_position` matrix.
    
    """
    
    # For every heavy atom id, get the column in the 
    # heavy_atom_position_matrix which holds its coords. Check that
    # these coords are the same as those held in the heavy confomer.    
    
    conf = mol.heavy_mol.GetConformer()
    pos_mat = mol.heavy_atom_position_matrix()    
    
    for atom_id in mol.heavy_ids:
        coord = np.array(conf.GetAtomPosition(atom_id))      
        
        column_i = mol.heavy_ids.index(atom_id)
        column = pos_mat.T[column_i]
        
        assert np.allclose(coord, column, atol=1e-8)    
    
def test_all_bonder_atom_distances():
    """
    Tests `all_bonder_atom_distances`.
    
    """

    # Check that the correct number of distances is found.
    assert sum(1 for _ in mol.all_heavy_atom_distances()) == 6
    
    s = sum(x for x, *_ in mol.all_heavy_atom_distances())
    assert round(s) == 33
    

def test_caching():
    
    # Use 2 .mol file to initialize 4 StructUnit instances such that
    # each .mol file initializes 2 StructUnit instances. Each pair 
    # initialzed from the same .mol file should be the same instance.
    # The ones from different .mol files should be different instnaces.
    
    bb_file = os.path.join('data', 'struct_unit', 
                           'struct_unit_caching_1.mol')
    lk_file = os.path.join('data', 'struct_unit', 
                           'struct_unit_caching_2.mol')

    bb = StructUnit(bb_file, minimal=True)
    bb2 = StructUnit(bb_file, minimal=True)
    
    lk = StructUnit(lk_file, minimal=True)
    lk2 = StructUnit(lk_file, minimal=True)
    
    assert bb is bb2
    assert bb is not lk
    assert bb is not lk2
    
    assert lk is lk2
    assert lk is not bb2
    assert lk2 is not bb2

def test_init():
    """
    Ensures that StructUnit instances are initiated correctly.
    
    This function only checks that the attributes exist and that
    their values of of the correct type. Checking that correct values
    are initialized is done in other tests.

    """     
    mol_file = os.path.join('data', 'struct_unit',
                            'struct_unit_init_aldehyde.mol')
    struct_unit = StructUnit(mol_file)
    
    assert isinstance(struct_unit.func_grp, FGInfo)
    assert isinstance(struct_unit.prist_mol, rdkit.Chem.rdchem.Mol)
    assert isinstance(struct_unit.heavy_mol, rdkit.Chem.rdchem.Mol)
    assert isinstance(struct_unit.prist_mol_file, str)
    assert isinstance(struct_unit.heavy_mol_file, str)
    assert isinstance(struct_unit.heavy_ids, list)

def test_make_atoms_heavy():
    """
    Tests `_make_atoms_heavy`.
    
    """
    
    mol_file = os.path.join('data', 'struct_unit',
                            'struct_unit_init_aldehyde.mol')
    struct_unit = StructUnit(mol_file)
    
    # Test that the position of the substituted atoms remains the same.
        
    i= 0
    for atom_id in flatten(struct_unit.find_functional_group_atoms()):
        atom = struct_unit.prist_mol.GetAtomWithIdx(atom_id)
        if (atom.GetAtomicNum() == 
            struct_unit.func_grp.target_atomic_num):
            prist_coord = struct_unit.atom_coords('prist', atom_id)
            heavy_coord = struct_unit.atom_coords('heavy',
                                               struct_unit.heavy_ids[i])
            i += 1
            assert np.array_equal(prist_coord, heavy_coord)  

def test_energy():
    """
    `energy` attribute must lazy and read the log file correctly.
    
    """
 
    # Load the StructUnit instance from the dump file.
    with open(obj_file_name, 'rb') as dump_file:
        struct_unit = pickle.load(dump_file)

    struct_unit.prist_mol_file = obj_file_name + '.mol'

    # `energy` should fail on unoptimized molecules as these will not 
    # have a .log file.
    with pytest.raises(AttributeError):
        struct_unit.energy

    # Check that the .log file is being read correctly.
    struct_unit.optimized = True
    assert struct_unit.energy == 1876.4573
    # Use another log file with a different energy to check that the
    # energy attribute is lazy. The log file getting changed but the
    # energy remaining the same proves this.
    struct_unit.prist_mol_file = obj_file_name + '2.mol'
    assert struct_unit.energy != -1876.4573
    # Delete the attribute so that the log file is reread and ensure
    # negative energies are gettting read correctly.
    delattr(struct_unit, 'energy')
    assert struct_unit.energy == -1876.4573    

def test_find_functional_group_atoms():
    """
    Make sure correct atoms are found in the functional groups.

    """
    
    # These are the expected atom ids for the test molecule.
    expected = ((1, 0, 12), (10, 11, 19))    
        
    func_grp_atoms = struct_unit.find_functional_group_atoms()

    assert func_grp_atoms == expected
 


def test_atom_distance_prist():
    """
    Test `atom_distance` when mol_type == 'prist'.
    
    """
    
    # Go through all combinations of atoms in the molecule. Calculate
    # the distance and compare it the distance calculated by the method.
    conf = struct_unit.prist_mol.GetConformer()
    for atom1, atom2 in it.combinations(
                                struct_unit.prist_mol.GetAtoms(), 2):
        atom1_id = atom1.GetIdx()
        atom2_id = atom2.GetIdx()
        assert (struct_unit.atom_distance('prist', 
                                          atom1_id, atom2_id) ==
               euclidean(conf.GetAtomPosition(atom1_id), 
                         conf.GetAtomPosition(atom2_id))) 
        
def test_atom_distance_heavy():
    """
    Test `atom_distance` when mol_type == 'heavy'.
    
    """
    
    # Go through all combinations of atoms in the molecule. Calculate
    # the distance and compare it the distance calculated by the method.
    conf = struct_unit.heavy_mol.GetConformer()
    for atom1, atom2 in it.combinations(
                                struct_unit.heavy_mol.GetAtoms(), 2):
        atom1_id = atom1.GetIdx()
        atom2_id = atom2.GetIdx()
        assert (struct_unit.atom_distance('heavy', 
                                          atom1_id, atom2_id) ==
               euclidean(conf.GetAtomPosition(atom1_id), 
                         conf.GetAtomPosition(atom2_id)))         

def test_centroid_functions_prist():
    """
    Tests functions related to centroid manipulation of prist molecule.
    
    Functions tested:
        > centroid
        > set_position
    
    """

    # Load the StructUnit instance from the dump file.
    with open(obj_file_name, 'rb') as dump_file:
        struct_unit = pickle.load(dump_file)
        
    # Get the centroid.
    prist_centroid = struct_unit.centroid('prist')
    # Get the heavy centroid to make sure its not changed in this test.    
    heavy_centroid = struct_unit.centroid('heavy')
    # Position the centroid.
    new_pos = np.array([10,15,25])
    struct_unit.set_position('prist', new_pos)
    # Check that the centroid is at the desired position and that it's
    # different to the original position.
    assert not np.allclose(prist_centroid, 
                           struct_unit.centroid('prist'), atol=1e-8)
    assert np.allclose(new_pos, struct_unit.centroid('prist'), 
                       atol = 1e-8)

    # Check that the heavy centroid is unmoved.
    assert np.array_equal(heavy_centroid, struct_unit.centroid('heavy'))
    assert not np.allclose(new_pos, struct_unit.centroid('heavy'), 
                       atol = 1e-8)

def test_centroid_functions_heavy():
    """
    Tests functions related to centroid manipulation of heavy molecule.
    
    Functions tested:
        > centroid
        > set_position
    
    """

    # Load the StructUnit instance from the dump file.
    with open(obj_file_name, 'rb') as dump_file:
        struct_unit = pickle.load(dump_file)
        
    # Get the centroid.
    heavy_centroid = struct_unit.centroid('heavy')
    # Get the prist centroid to make sure its not changed in this test.    
    prist_centroid = struct_unit.centroid('prist')
    # Position the centroid.
    new_pos = np.array([10,15,25])
    struct_unit.set_position('heavy', new_pos)
    # Check that the centroid is at the desired position and that it's
    # different to the original position.
    assert not np.allclose(heavy_centroid, 
                           struct_unit.centroid('heavy'), atol=1e-8)
    assert np.allclose(new_pos, struct_unit.centroid('heavy'), 
                       atol = 1e-8)

    # Check that the prist centroid is unmoved.
    assert np.array_equal(prist_centroid, struct_unit.centroid('prist'))
    assert not np.allclose(new_pos, struct_unit.centroid('prist'), 
                       atol = 1e-8) 
 
def test_shift():
    """
    Test `shift`.      
    
    """
    
    # Load the StructUnit instance from the dump file.
    with open(obj_file_name, 'rb') as dump_file:
        struct_unit = pickle.load(dump_file)

    # Shifting the same StructUnit twice should return two 
    # ``rdkit.Chem.rdchem.Mol`` instances with conformers describing the
    # same atomic positions. Furthermore, the original conformer in the
    # ``StructUnit`` instance should be unchanged.
    og_conformer = struct_unit.prist_mol.GetConformer()

    shift = [10,10,10]    
    a = struct_unit.shift('prist', shift)
    a_conformer = a.GetConformer()
    
    b = struct_unit.shift('prist', shift)
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
        assert atom1_coords.x == og_coords.x + shift[0]
        assert atom1_coords.y == og_coords.y + shift[0]
        assert atom1_coords.z == og_coords.z + shift[0]
    
def test_set_position_from_matrix_prist():
    """
    Tests `set_position_from_matrix` when mol_type == 'prist'.
    
    """
    
    # Make a position matrix where each coordinate is set to [1,2,3] for
    # every atom. Set this as the position matrix for the molecule. Make
    # sure that only the `prist` molecule is modified.
    
    # Load the StructUnit instance from the dump file.
    with open(obj_file_name, 'rb') as dump_file:
        struct_unit = pickle.load(dump_file) 
    
    pos_mat = []
    for x in range(struct_unit.prist_mol.GetNumAtoms()):
        pos_mat.append([1,2,3])
    
    pos_mat = np.matrix(pos_mat).T

    struct_unit.set_position_from_matrix('prist', pos_mat)
    
    for _, coord in struct_unit.all_atom_coords('prist'):
        assert np.array_equal(coord, [1,2,3])

    for _, coord in struct_unit.all_atom_coords('heavy'):
        assert not np.array_equal(coord, [1,2,3])
    
    # Centroid should also be in [1,2,3].
    assert np.array_equal(struct_unit.centroid('prist'), [1,2,3])    

def test_set_position_from_matrix_heavy():
    """
    Tests `set_position_from_matrix` when mol_type == 'heavy'.
    
    """
    
    # Make a position matrix where each coordinate is set to [1,2,3] for
    # every atom. Set this as the position matrix for the molecule. Make
    # sure that only the `heavy` molecule is modified.
    
    # Load the StructUnit instance from the dump file.
    with open(obj_file_name, 'rb') as dump_file:
        struct_unit = pickle.load(dump_file) 
    
    pos_mat = []
    for x in range(struct_unit.heavy_mol.GetNumAtoms()):
        pos_mat.append([1,2,3])
    
    pos_mat = np.matrix(pos_mat).T

    struct_unit.set_position_from_matrix('heavy', pos_mat)
    
    for _, coord in struct_unit.all_atom_coords('heavy'):
        assert np.array_equal(coord, [1,2,3])

    for _, coord in struct_unit.all_atom_coords('prist'):
        assert not np.array_equal(coord, [1,2,3])
    
    # Centroid should also be in [1,2,3].
    assert np.array_equal(struct_unit.centroid('heavy'), [1,2,3])

def test_set_heavy_mol_orientation():
    """
    Tests `_set_heavy_mol_orientation`.
    
    """

    # Takes the first direction vector between two heavy atoms and
    # sets it to [1,2,3]. Ensure that only heavy molecule is modified
    # and that all coordinates are shifted.

    # Load the StructUnit instance from the dump file.
    with open(obj_file_name, 'rb') as dump_file:
        struct_unit = pickle.load(dump_file) 
    
    # Save the position matrix for comparison later.
    prist_pos_mat = struct_unit.position_matrix('prist')
    
    # Set orienation.
    struct_unit._set_heavy_mol_orientation(
    next(struct_unit.heavy_direction_vectors()), [1,2,3])
    
    # Check correct orientation.
    assert np.allclose(next(struct_unit.heavy_direction_vectors()), 
                       normalize_vector([1,2,3]),
                       atol=1e-3) 
    # Check that pristine molecule was not modified.
    assert np.array_equal(prist_pos_mat, 
                          struct_unit.position_matrix('prist'))



        


def test_centroid_centroid_dir_vector():
    """
    Tests `centroid_centroid_dir_vector`.
    
    """
    
    mol_centroid = struct_unit.centroid('heavy')
    h_atom_centroid =struct_unit.heavy_atom_centroid()
    
    assert np.array_equal(
        normalize_vector(mol_centroid - h_atom_centroid),
                          struct_unit.centroid_centroid_dir_vector())

def test_heavy_atom_centroid_functions():
    """
    Tests functions which manipulate the ``heavy_atom_centroid``.
    
    Testsed functions include:
         > heavy_atom_centroid
         > set_heavy_atom_centroid
    
    """
    
    # Load the StructUnit instance from the dump file.
    with open(obj_file_name, 'rb') as dump_file:
        struct_unit = pickle.load(dump_file)     
    
    # Get the heavy atom centroid. Then move it. Check that the prist
    # molecule was unchanged and that the move was successful.    
    
    prist_centroid = struct_unit.centroid('prist')
    heavy_centroid = struct_unit.centroid('heavy')
    heavy_atom_centroid = struct_unit.heavy_atom_centroid()    
    
    # Move and check.    
    struct_unit.set_heavy_atom_centroid([1,2,3])
    assert np.allclose(struct_unit.heavy_atom_centroid(), 
                       [1,2,3], atol=1e-8)
    assert np.array_equal(prist_centroid, struct_unit.centroid('prist'))
    assert not np.allclose(heavy_centroid, 
                           struct_unit.centroid('heavy'), atol=1e-4)
    assert not np.allclose(heavy_atom_centroid,
                           struct_unit.heavy_atom_centroid(), atol=1e-4)

def test_write_mol_file():
    """
    Tests `write_mol_file`.
    
    """
    
    # Load the StructUnit instance from the dump file.
    with open(obj_file_name, 'rb') as dump_file:
        struct_unit = pickle.load(dump_file)    
    
    struct_unit.prist_mol_file = 'delete_this.mol'
    struct_unit.heavy_mol_file = 'delete_this_heavy.mol' 
    
    struct_unit.write_mol_file('prist')
    struct_unit.write_mol_file('heavy')
    
    
    # Get the expected output as a string.
    prist_name = os.path.join('data','struct_unit', 
                              'write_test_prist.mol')
    heavy_name = os.path.join('data','struct_unit', 
                              'write_test_heavy.mol')  
    
    with open(prist_name, 'r') as prist_file:
        exp_output_prist = prist_file.read()
    
    with open(heavy_name, 'r') as heavy_file:
        exp_output_heavy = heavy_file.read() 
    
    # Get the written output as a string.
    with open('delete_this.mol', 'r') as out_file:
        output_prist = out_file.read()
        
    with open('delete_this_heavy.mol', 'r') as out_file:
        output_heavy = out_file.read()
        
    assert exp_output_prist == output_prist
    assert exp_output_heavy == output_heavy
    
    
    