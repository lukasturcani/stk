import os
import pickle
import pytest
import numpy as np
import itertools as it
from scipy.spatial.distance import euclidean

from ...classes import (MacroMolecule, FourPlusSix, StructUnit3, 
                        StructUnit2, EightPlusTwelve)

from ...convenience_functions import periodic_table

obj_file_name = os.path.join('data', 'macromolecule',
                             'macromolecule_test_obj')

bb_file = os.path.join('data', 'macromolecule', 
                       'macromolecule_bb_amine.mol')
lk_file = os.path.join('data', 'macromolecule', 
                       'macromolecule_lk_aldehyde_1.mol')

bb = StructUnit3(bb_file)
lk = StructUnit2(lk_file)    
building_blocks = (bb, lk)
mol = MacroMolecule(building_blocks, FourPlusSix,
                    'delete_this_1.mol')

def test_caching():
    """
    Cages created with same arguments should return the same instance.

    """
    
    # First init a new MacroMolecule using the same values as prviously.
    # This should yield the same MacroMolecule.
    mol2 = MacroMolecule(building_blocks, FourPlusSix, 
                                 'delete_this_2.mol')    
    assert mol is mol2

    # Change the order of the building blocks in the tuple, this
    # should still be the same MacroMolecule.
    building_blocks2 = (lk, bb)
    mol2 = MacroMolecule(building_blocks2, FourPlusSix, 
                         'delete_this_3.mol')    
    assert mol is  mol2
    
    # Change the topology, this should make a new MacroMolecule.     
    mol2 = MacroMolecule(building_blocks, EightPlusTwelve,
                         'delete_this_4.mol')    
    assert mol is not mol2
    
    
    # Change one of the building blocks, this should make a new 
    # MacroMolecule
    lk2_file = os.path.join('data', 'macromolecule',
                            'macromolecule_lk_aldehyde_2.mol')
                                        
    lk2 = StructUnit2(lk2_file)
    building_blocks2 = (bb, lk2)
    mol2 = MacroMolecule(building_blocks2, FourPlusSix, 
                          'delete_this_5.mol')

    assert mol is not mol2

def test_energy():
    """
    `energy` attribute must lazy and read the log file correctly.
    
    """
 
    # Load the StructUnit instance from the dump file.
    with open(obj_file_name, 'rb') as dump_file:
        macro_mol = pickle.load(dump_file)

    macro_mol.prist_mol_file = obj_file_name + '.mol'

    # `energy` should fail on unoptimized molecules as these will not 
    # have a .log file.
    with pytest.raises(AttributeError):
        macro_mol.energy

    # Check that the .log file is being read correctly.
    macro_mol.optimized = True
    assert macro_mol.energy == 1876.4573
    # Use another log file with a different energy to check that the
    # energy attribute is lazy. The log file getting changed but the
    # energy remaining the same proves this.
    macro_mol.prist_mol_file = obj_file_name + '2.mol'
    assert macro_mol.energy != -1876.4573
    # Delete the attribute so that the log file is reread and ensure
    # negative energies are gettting read correctly.
    delattr(macro_mol, 'energy')
    assert macro_mol.energy == -1876.4573  

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
    prist_name = os.path.join('data','macromolecule', 
                              'write_test_prist.mol')
    heavy_name = os.path.join('data','macromolecule', 
                              'write_test_heavy.mol')  
                
    from itertools import zip_longest
    with open(prist_name, 'r') as prist_file:
        
        with open(struct_unit.prist_mol_file, 'r') as out_file:
            for l1, l2 in zip_longest(prist_file, out_file):
                assert l1 == l2
   
   
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

def test_graph():
    """
    Tests the output of the `graph` method.    
    
    """
    # Test the pristine version first.
    graph = mol.graph('prist')
    expected_nodes = 236
    expected_edges = 256
    assert len(graph.nodes()) == expected_nodes
    assert len(graph.edges()) == expected_edges

    
    # Test the heavy version second.
    graph = mol.graph('heavy')
    expected_nodes = 224
    expected_edges = 244

    assert len(graph.nodes()) == expected_nodes
    assert len(graph.edges()) == expected_edges

def test_atom_symbol():
    """
    Tests the `atom_symbol` method.
    
    """
    
    # First test when `mol_type` == 'prist'.
    for atom in mol.prist_mol.GetAtoms():
        atom_id = atom.GetIdx()
        atom_sym = periodic_table[atom.GetAtomicNum()]
        assert atom_sym == mol.atom_symbol('prist', atom_id)
    
    # Test when `mol_typ` == 'heavy'.
    for atom in mol.heavy_mol.GetAtoms():
        atom_id = atom.GetIdx()
        atom_sym = periodic_table[atom.GetAtomicNum()]
        assert atom_sym == mol.atom_symbol('heavy', atom_id)
        
def test_atom_coords_prist():
    """
    Tests `atom_coords` when mol_type == 'prist'.
    
    """
    
    conf = mol.prist_mol.GetConformer()    
    
    for atom in mol.prist_mol.GetAtoms():
        atom_id = atom.GetIdx()
        coords = mol.atom_coords('prist', atom_id)
        conf_coords = conf.GetAtomPosition(atom_id)
        assert np.allclose(coords, conf_coords, atol=1e-8)
    
def test_atom_coords_heavy():
    """
    Test `atom_coords` when mol_type == 'heavy'.
    
    """
    
    conf = mol.heavy_mol.GetConformer()    
    
    for atom in mol.heavy_mol.GetAtoms():
        atom_id = atom.GetIdx()
        coords = mol.atom_coords('heavy', atom_id)
        conf_coords = conf.GetAtomPosition(atom_id)
        assert np.allclose(coords, conf_coords, atol=1e-8) 

def test_all_atom_coords_prist():
    """
    Test `all_atom_coords` when mol_type == 'prist'.

    """
        
    conf = mol.prist_mol.GetConformer()
    for (atom_id, coord), atom in it.zip_longest(
                                mol.all_atom_coords('prist'), 
                                mol.prist_mol.GetAtoms()):
        
        assert atom_id == atom.GetIdx()
        conf_coord = np.array(conf.GetAtomPosition(atom_id))
        assert np.allclose(coord, conf_coord, atol=1e-8)
        
def test_all_atom_coords_heavy():
    """
    Test `all_atom_coords` when mol_type == 'heavy'.

    """
    
    # This test takes too long if all distances are checked. Randomly 
    # sample 10 instead.    
    
    conf = mol.heavy_mol.GetConformer()
    for (atom_id, coord), atom in it.zip_longest(
                                mol.all_atom_coords('heavy'), 
                                mol.heavy_mol.GetAtoms()):
        
        assert atom_id == atom.GetIdx()
        conf_coord = np.array(conf.GetAtomPosition(atom_id))
        assert np.allclose(coord, conf_coord, atol=1e-8)            

def test_atom_distance_prist():
    """
    Test `atom_distance` when mol_type == 'prist'.
    
    """

    # This test takes too long if all distances are checked. Randomly 
    # sample 10 instead.  
    
    # Go through all combinations of atoms in the molecule. Calculate
    # the distance and compare it the distance calculated by the method.
    conf = mol.prist_mol.GetConformer()
    for atom1, atom2 in it.combinations(
                          np.random.choice(
                                mol.prist_mol.GetAtoms(),size=10), 2):
        atom1_id = atom1.GetIdx()
        atom2_id = atom2.GetIdx()
        assert (mol.atom_distance('prist', 
                                          atom1_id, atom2_id) ==
               euclidean(conf.GetAtomPosition(atom1_id), 
                         conf.GetAtomPosition(atom2_id))) 
        
def test_atom_distance_heavy():
    """
    Test `atom_distance` when mol_type == 'heavy'.
    
    """

    # This test takes too long if all distances are checked. Randomly 
    # sample 10 instead. 
    
    # Go through all combinations of atoms in the molecule. Calculate
    # the distance and compare it the distance calculated by the method.
    conf = mol.heavy_mol.GetConformer()
    for atom1, atom2 in it.combinations(
                          np.random.choice(
                              mol.heavy_mol.GetAtoms(), size=10), 2):
        atom1_id = atom1.GetIdx()
        atom2_id = atom2.GetIdx()
        assert (mol.atom_distance('heavy', 
                                          atom1_id, atom2_id) ==
               euclidean(conf.GetAtomPosition(atom1_id), 
                         conf.GetAtomPosition(atom2_id)))   
   
def test_all_heavy_atom_distances():
    """
    Tests `all_heavy_atom_distances`.
    
    """

    # Check that the correct number of distances is found.
    assert sum(1 for _ in mol.all_heavy_atom_distances()) == 276

def test_centroid_functions_prist():
    """
    Tests functions related to centroid manipulation of prist molecule.
    
    Functions tested:
        > centroid
        > set_position
    
    """

    # Load the StructUnit instance from the dump file.
    with open(obj_file_name, 'rb') as dump_file:
        mol = pickle.load(dump_file)
        
    # Get the centroid.
    prist_centroid = mol.centroid('prist')
    # Get the heavy centroid to make sure its not changed in this test.    
    heavy_centroid = mol.centroid('heavy')
    # Position the centroid.
    new_pos = np.array([10,15,25])
    mol.set_position('prist', new_pos)
    # Check that the centroid is at the desired position and that it's
    # different to the original position.
    assert not np.allclose(prist_centroid, 
                           mol.centroid('prist'), atol=1e-8)
    assert np.allclose(new_pos, mol.centroid('prist'), 
                       atol = 1e-8)

    # Check that the heavy centroid is unmoved.
    assert np.array_equal(heavy_centroid, mol.centroid('heavy'))
    assert not np.allclose(new_pos, mol.centroid('heavy'), 
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
        mol = pickle.load(dump_file)
        
    # Get the centroid.
    heavy_centroid = mol.centroid('heavy')
    # Get the prist centroid to make sure its not changed in this test.    
    prist_centroid = mol.centroid('prist')
    # Position the centroid.
    new_pos = np.array([10,15,25])
    mol.set_position('heavy', new_pos)
    # Check that the centroid is at the desired position and that it's
    # different to the original position.
    assert not np.allclose(heavy_centroid, 
                           mol.centroid('heavy'), atol=1e-8)
    assert np.allclose(new_pos, mol.centroid('heavy'), 
                       atol = 1e-8)

    # Check that the prist centroid is unmoved.
    assert np.array_equal(prist_centroid, mol.centroid('prist'))
    assert not np.allclose(new_pos, mol.centroid('prist'), 
                       atol = 1e-8) 

def test_center_of_mass():
    """
    Tests `center_of_mass`.
    
    """
    
    # Calculate the center of mass.
    coord_sum = 0
    total_mass = 0
    for atom_id, coord in mol.all_atom_coords('prist'):
        atom_mass = mol.prist_mol.GetAtomWithIdx(atom_id).GetMass()
        total_mass += atom_mass
        scaled_coord = np.multiply(atom_mass, coord)
        coord_sum = np.add(scaled_coord, coord_sum)

    com = np.divide(coord_sum, total_mass)
    assert np.allclose(mol.center_of_mass('prist'), com, atol=1e-6)

    coord_sum = 0
    total_mass = 0
    for atom_id, coord in mol.all_atom_coords('heavy'):
        atom_mass = mol.heavy_mol.GetAtomWithIdx(atom_id).GetMass()
        total_mass += atom_mass
        scaled_coord = np.multiply(atom_mass, coord)
        coord_sum = np.add(scaled_coord, coord_sum)

    com = np.divide(coord_sum, total_mass)
    assert np.allclose(mol.center_of_mass('heavy'), com, atol=1e-6)

def test_shift():
    """
    Test `shift`.      
    
    """
    
    # Load the StructUnit instance from the dump file.
    with open(obj_file_name, 'rb') as dump_file:
        macro_mol = pickle.load(dump_file)

    # Shifting the same StructUnit twice should return two 
    # ``rdkit.Chem.rdchem.Mol`` instances with conformers describing the
    # same atomic positions. Furthermore, the original conformer in the
    # ``StructUnit`` instance should be unchanged.
    og_conformer = macro_mol.prist_mol.GetConformer()

    shift = [10,10,10]    
    a = macro_mol.shift('prist', shift)
    a_conformer = a.GetConformer()
    
    b = macro_mol.shift('prist', shift)
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

def test_same():
    """
    Tests the `same_cage` method.    
    
    Cages initialized from the same arguments should return ``True`` 
    through this method, even if the ``Cage`` class stops being cached.    
    
    """
    
    a = MacroMolecule.testing_init('a', 'b', 'c')
    b = MacroMolecule.testing_init('a', 'a', 'b')
    c = MacroMolecule.testing_init('a', 'a', 'b')
    d = MacroMolecule.testing_init('a', 'b', 'a')

    assert not a.same(b)
    assert b.same(c)    
    assert c.same(b)
    assert not d.same(c)
        
def test_comparison():
    """
    Checks ``==``, ``>``, ``>=``, etc. operators.    
    
    """
    
    # Generate cages with various fitnesses.
    a = MacroMolecule.testing_init('a','a','a')
    a.fitness = 1
    
    b = MacroMolecule.testing_init('b', 'b', 'b')
    b.fitness = 1
    
    c = MacroMolecule.testing_init('c', 'c', 'c')
    c.fitness = 2
    
    # Comparison operators should compare their fitness.
    assert not a < b
    assert a <= b
    assert a == b
    assert c > b
    assert c >= a
