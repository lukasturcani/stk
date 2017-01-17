import rdkit.Chem as chem
from os.path import join
import itertools as it
import numpy as np
from scipy.spatial.distance import euclidean


from ..classes import Molecule
from ..convenience_tools import periodic_table

data_dir = join('data', 'molecule')

mol = Molecule.__new__(Molecule)
mol.prist_mol = chem.MolFromMolFile(join(data_dir, 'prist.mol'),
                                    removeHs=False, sanitize=False)
mol.heavy_mol = chem.MolFromMolFile(join(data_dir, 'heavy.mol'),
                                    removeHs=False, sanitize=False)
mol.heavy_ids = [0, 9, 12, 15]

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
        
    conf = mol.heavy_mol.GetConformer()
    for (atom_id, coord), atom in it.zip_longest(
                                mol.all_atom_coords('heavy'), 
                                mol.heavy_mol.GetAtoms()):
        
        assert atom_id == atom.GetIdx()
        conf_coord = np.array(conf.GetAtomPosition(atom_id))
        assert np.allclose(coord, conf_coord, atol=1e-8)
        
def test_all_heavy_atom_distances():
    """
    Tests `all_heavy_atom_distances`.
    
    """

    # Check that the correct number of distances is found.
    assert sum(1 for _ in mol.all_heavy_atom_distances()) == 6
    
    s = sum(x for x, *_ in mol.all_heavy_atom_distances())
    assert round(s) == 33
    
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
    
def test_atom_distance_prist():
    """
    Test `atom_distance` when mol_type == 'prist'.
    
    """
    
    # Go through all combinations of atoms in the molecule. Calculate
    # the distance and compare it the distance calculated by the method.
    conf = mol.prist_mol.GetConformer()
    for atom1, atom2 in it.combinations(mol.prist_mol.GetAtoms(), 2):
        atom1_id = atom1.GetIdx()
        atom2_id = atom2.GetIdx()
        assert (mol.atom_distance('prist', atom1_id, atom2_id) ==
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
    for atom1, atom2 in it.combinations(mol.heavy_mol.GetAtoms(), 2):
        atom1_id = atom1.GetIdx()
        atom2_id = atom2.GetIdx()
        assert (mol.atom_distance('heavy', atom1_id, atom2_id) ==
                euclidean(conf.GetAtomPosition(atom1_id), 
                          conf.GetAtomPosition(atom2_id))) 
    
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
    
def test_centroid_functions_prist():
    """
    Tests functions related to centroid manipulation of prist molecule.
    
    Functions tested:
        > centroid
        > set_position
    
    """
        
    # Get the centroid.
    prist_centroid = mol.centroid('prist')
    # Get the heavy centroid to make sure its not changed in this test.    
    heavy_centroid = mol.centroid('heavy')
    # Position the centroid.
    new_pos = np.array([25,15,10])
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
        
    # Get the centroid.
    heavy_centroid = mol.centroid('heavy')
    # Get the prist centroid to make sure its not changed in this test.    
    prist_centroid = mol.centroid('prist')
    # Position the centroid.
    new_pos = np.array([10.,15.,25.])
    mol.set_position('heavy', new_pos)
    # Check that the centroid is at the desired position and that it's
    # different to the original position.
    assert not np.allclose(heavy_centroid, 
                           mol.centroid('heavy'), atol=1e-8)
    assert np.allclose(new_pos, mol.centroid('heavy'), 
                       atol = 1e-8)

    # Check that the prist centroid is unmoved.
    assert np.allclose(prist_centroid, mol.centroid('prist'), atol=1e-8)
    assert not np.allclose(new_pos, mol.centroid('prist'), atol = 1e-8)
    
def test_graph():
    """
    Tests the output of the `graph` method.    
    
    """
    # Test the pristine version first.
    graph = mol.graph('prist')
    expected_nodes = 24
    expected_edges = 24
    assert len(graph.nodes()) == expected_nodes
    assert len(graph.edges()) == expected_edges

    
    # Test the heavy version second.
    graph = mol.graph('heavy')
    expected_nodes = 16
    expected_edges = 16

    assert len(graph.nodes()) == expected_nodes
    assert len(graph.edges()) == expected_edges
    
def test_heavy_atom_position_matrix():
    """
    Tests the output of `heavy_atom_position` matrix.
    
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
        
def test_position_matrix_prist():
    """
    Test `postion_matrix` when mol_type == 'prist'.
    
    """
    
    # Go through each atom id. For each atom id get the column in the 
    # position matrix with that id as its index. Make sure that the data
    # is the same. 
    pos_mat1 = mol.position_matrix('prist')
    conf = mol.prist_mol.GetConformer()
       
    for atom in mol.prist_mol.GetAtoms():
        atom_id = atom.GetIdx()
        conf_coord = np.array(conf.GetAtomPosition(atom_id))   
        mat_coord = pos_mat1.T[atom_id]

        assert np.allclose(conf_coord, mat_coord, atol = 1e-8)

def test_position_matrix_heavy():
    """
    Test `postion_matrix` when mol_type == 'heavy'.
    
    """
    
    # Go through each atom id. For each atom id get the column in the 
    # position matrix with that id as its index. Make sure that the data
    # is the same. 
    pos_mat1 = mol.position_matrix('heavy')
    conf = mol.heavy_mol.GetConformer()
       
    for atom in mol.heavy_mol.GetAtoms():
        atom_id = atom.GetIdx()
        conf_coord = np.array(conf.GetAtomPosition(atom_id))   
        mat_coord = pos_mat1.T[atom_id]

        assert np.allclose(conf_coord, mat_coord, atol = 1e-8)

def test_rotate():
    prist_pos_mat = np.matrix([[2.492382132029387520e+01, 2.529570505995089036e+01, 2.530589753727284474e+01, 2.561980255444483845e+01, 2.703139356920017278e+01, 2.470134676063525347e+01, 2.468318917529596135e+01, 2.372283110870242950e+01, 2.433590679503502585e+01, 2.398086135664957652e+01, 2.489387059221622422e+01, 2.630704756049889426e+01, 2.432154250041805454e+01, 2.602499611683657932e+01, 2.543308750983269917e+01, 2.761824314389179591e+01, 2.729016888784238404e+01, 2.499823142630443229e+01, 2.368347697672986030e+01, 2.567282839046659504e+01, 2.278356608270766515e+01, 2.368851529192009409e+01, 2.330210749325866715e+01, 2.438156278959523249e+01],
                               [1.768830014659173244e+01, 1.631106674498808928e+01, 1.545619308055846552e+01, 1.397723849756132886e+01, 1.377422828201651406e+01, 1.341339980610046645e+01, 1.427136558384375853e+01, 1.371962211690825306e+01, 1.573139958145930706e+01, 1.771337616367322454e+01, 1.825369393352212555e+01, 1.633978675506987699e+01, 1.551240081079739497e+01, 1.586411641831839248e+01, 1.341159534980574897e+01, 1.415239038290292051e+01, 1.427845771863169588e+01, 1.238450683382546913e+01, 1.334688177349772609e+01, 1.423696999519013318e+01, 1.372972150174510730e+01, 1.430470755596099330e+01, 1.579540944679532544e+01, 1.633317152023594332e+01],
                               [1.052050853501445538e+01, 1.018688404302633010e+01, 1.146046174864169487e+01, 1.118633737114347504e+01, 1.085072680807397916e+01, 1.009054230807505093e+01, 8.818624920661767064e+00, 7.858930342041492878e+00, 9.139953011396050542e+00, 1.090799298474652446e+01, 9.672311736722228659e+00, 9.764578837643462350e+00, 1.194652185402330069e+01, 1.218328681237754196e+01, 1.210718497441264319e+01, 1.159473913376732135e+01, 1.000543072056803062e+01, 9.846837121360788814e+00, 1.050015379090626233e+01, 8.347829329597315606e+00, 8.254533350678615378e+00, 7.025094638550643644e+00, 9.507896149462199631e+00, 8.222639477108881323e+00]])
    
    mol.rotate('prist', np.pi/3.2, [1,2,3])
    assert np.allclose(mol.position_matrix('prist'), prist_pos_mat, atol=1e-8)
    
