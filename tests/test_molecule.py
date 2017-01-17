import rdkit.Chem as chem
from os.path import join
import itertools as it
import numpy as np
from scipy.spatial.distance import euclidean


from ..classes import Molecule
from ..convenience_tools import periodic_table

data_dir = join('data', 'molecule')

mol = Molecule.__new__(Molecule)
mol.mol = chem.MolFromMolFile(join(data_dir, 'molecule.mol'),
                                    removeHs=False, sanitize=False)
og = mol.position_matrix()
def test_all_atom_coords():
    """
    Test `all_atom_coords`.

    """
        
    conf = mol.mol.GetConformer()
    for (atom_id, coord), atom in it.zip_longest(
                                mol.all_atom_coords(), 
                                mol.mol.GetAtoms()):
        
        assert atom_id == atom.GetIdx()
        conf_coord = np.array(conf.GetAtomPosition(atom_id))
        assert np.allclose(coord, conf_coord, atol=1e-8)
        
def test_atom_coords():
    """
    Tests `atom_coords`.
    
    """
    
    conf = mol.mol.GetConformer()    
    
    for atom in mol.mol.GetAtoms():
        atom_id = atom.GetIdx()
        coords = mol.atom_coords(atom_id)
        conf_coords = conf.GetAtomPosition(atom_id)
        assert np.allclose(coords, conf_coords, atol=1e-8)  
    
def test_atom_distance():
    """
    Test `atom_distance`.
    
    """
    
    # Go through all combinations of atoms in the molecule. Calculate
    # the distance and compare it the distance calculated by the method.
    conf = mol.mol.GetConformer()
    for atom1, atom2 in it.combinations(mol.mol.GetAtoms(), 2):
        atom1_id = atom1.GetIdx()
        atom2_id = atom2.GetIdx()
        assert (mol.atom_distance(atom1_id, atom2_id) ==
               euclidean(conf.GetAtomPosition(atom1_id), 
                         conf.GetAtomPosition(atom2_id))) 
    
def test_atom_symbol():
    """
    Tests the `atom_symbol` method.
    
    """
    
    for atom in mol.mol.GetAtoms():
        atom_id = atom.GetIdx()
        atom_sym = periodic_table[atom.GetAtomicNum()]
        assert atom_sym == mol.atom_symbol(atom_id)

        
def test_center_of_mass():
    """
    Tests `center_of_mass`.
    
    """
    
    # Calculate the center of mass.
    coord_sum = 0
    total_mass = 0
    for atom_id, coord in mol.all_atom_coords():
        atom_mass = mol.mol.GetAtomWithIdx(atom_id).GetMass()
        total_mass += atom_mass
        scaled_coord = np.multiply(atom_mass, coord)
        coord_sum = np.add(scaled_coord, coord_sum)

    com = np.divide(coord_sum, total_mass)
    assert np.allclose(mol.center_of_mass(), com, atol=1e-6)
    
def test_centroid_functions():
    """
    Tests functions related to centroid manipulation of the molecule.
    
    Functions tested:
        > centroid
        > set_position
    
    """
        
    try:
        # Get the centroid.
        prist_centroid = mol.centroid()
        # Position the centroid.
        new_pos = np.array([25,15,10])
        mol.set_position(new_pos)
        # Check that the centroid is at the desired position and that it's
        # different to the original position.
        assert not np.allclose(prist_centroid, 
                               mol.centroid(), atol=1e-8)
        assert np.allclose(new_pos, mol.centroid(), 
                           atol = 1e-8)
                           
    except Exception as ex:
        mol.set_position_from_matrix(og)
        raise ex

def test_graph():
    """
    Tests the output of the `graph` method.    
    
    """
    # Test the pristine version first.
    graph = mol.graph()
    expected_nodes = 24
    expected_edges = 24
    assert len(graph.nodes()) == expected_nodes
    assert len(graph.edges()) == expected_edges
        
def test_position_matrix():
    """
    Test `postion_matrix`.
    
    """
    
    # Go through each atom id. For each atom id get the column in the 
    # position matrix with that id as its index. Make sure that the data
    # is the same. 
    pos_mat1 = mol.position_matrix()
    conf = mol.mol.GetConformer()
       
    for atom in mol.mol.GetAtoms():
        atom_id = atom.GetIdx()
        conf_coord = np.array(conf.GetAtomPosition(atom_id))   
        mat_coord = pos_mat1.T[atom_id]

        assert np.allclose(conf_coord, mat_coord, atol = 1e-8)

def test_set_position_from_matrix():
    try:
        new_pos_mat = np.matrix([[0 for x in range(3)] for y in 
                                        range(mol.mol.GetNumAtoms())])
        mol.set_position_from_matrix(new_pos_mat.T)
        for _, atom_coord in mol.all_atom_coords():
            assert np.allclose(atom_coord, [0,0,0], atol=1e-8) 
    except Exception as ex:
        mol.set_position_from_matrix(og)
        raise ex

def test_shift():
    s= np.array([10,-20,5])
    mol2 = mol.shift(s)
    conf = mol2.GetConformer()
    for atom in mol2.GetAtoms():
        atomid = atom.GetIdx()
        pos = conf.GetAtomPosition(atomid)
        should_be = mol.atom_coords(atomid) + s
        assert np.allclose(should_be, pos,atol=1e-8)