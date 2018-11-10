from os.path import join
import itertools as it
import numpy as np
from scipy.spatial.distance import euclidean
import stk


def test_all_atom_coords(mol):
    """
    Test `all_atom_coords`.

    """

    natoms = mol.mol.GetNumAtoms()
    for i, (atom, coord) in enumerate(mol.all_atom_coords(), 1):
        assert atom < natoms
        assert type(atom) == int
        assert type(coord) == np.ndarray
        assert len(coord) == 3
    assert natoms == i


def test_atom_coords(mol):
    """
    Tests `atom_coords`.

    """

    conf = mol.mol.GetConformer()

    for atom in mol.mol.GetAtoms():
        atom_id = atom.GetIdx()
        coords = mol.atom_coords(atom_id)
        conf_coords = conf.GetAtomPosition(atom_id)
        assert np.allclose(coords, conf_coords, atol=1e-8)


def test_atom_distance(mol):
    """
    Test `atom_distance`.

    """

    # Go through all combinations of atoms in the molecule. Calculate
    # the distance and compare it the distance calculated by the
    # method.
    conf = mol.mol.GetConformer()
    for atom1, atom2 in it.combinations(mol.mol.GetAtoms(), 2):
        atom1_id = atom1.GetIdx()
        atom2_id = atom2.GetIdx()
        assert (mol.atom_distance(atom1_id, atom2_id) ==
                euclidean(conf.GetAtomPosition(atom1_id),
                          conf.GetAtomPosition(atom2_id)))


def test_atom_symbol(mol):
    """
    Tests the `atom_symbol` method.

    """

    for atom in mol.mol.GetAtoms():
        atom_id = atom.GetIdx()
        atom_sym = stk.periodic_table[atom.GetAtomicNum()]
        assert atom_sym == mol.atom_symbol(atom_id)


def test_cavity_size(cc3):
    assert abs(cc3.cavity_size()-6.30569486) < 1e-4


def test_center_of_mass(mol):
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


def test_centroid_functions(tmp_mol):
    """
    Tests functions related to centroid manipulation of the molecule.

    Functions tested:
        > centroid
        > set_position

    """

    # Save the coordinates of the new centroid.
    new_centroid = tmp_mol.centroid() + np.array([10, 20, 4])
    tmp_mol.set_position(new_centroid)
    # Check that the centroid is at the desired position.
    assert np.allclose(new_centroid, tmp_mol.centroid(), atol=1e-8)


def test_graph(mol):
    """
    Tests the output of the `graph` method.

    """

    graph = mol.graph()
    assert len(graph.nodes()) == mol.mol.GetNumAtoms()
    assert len(graph.edges()) == mol.mol.GetNumBonds()


def test_max_diameter(tmp_mol):
        # Make a position matrix which sets all atoms to the origin except
        # 2 and 13. These should be placed a distance of 100 apart.
        pos_mat = [[0 for x in range(3)] for
                   y in range(tmp_mol.mol.GetNumAtoms())]
        pos_mat[1] = [0, -50, 0]
        pos_mat[12] = [0, 50, 0]
        tmp_mol.set_position_from_matrix(np.matrix(pos_mat).T)

        d, id1, id2 = tmp_mol.max_diameter()
        # Note that it is not exactly 100 because of the Van der Waals
        # radii of the atoms.
        assert d > 100 and d < 105
        assert id1 == 1
        assert id2 == 12


def test_position_matrix(mol):
    """
    Test `postion_matrix`.

    """

    # Go through each atom id. For each atom id get the column in the
    # position matrix with that id as its index. Make sure that the
    # data is the same.
    pos_mat1 = mol.position_matrix()
    conf = mol.mol.GetConformer()

    for atom in mol.mol.GetAtoms():
        atom_id = atom.GetIdx()
        conf_coord = np.array(conf.GetAtomPosition(atom_id))
        mat_coord = pos_mat1.T[atom_id]

        assert np.allclose(conf_coord, mat_coord, atol=1e-8)


def test_same(mol, tmp_mol, mol2):
    """
    Tests the `same()` method.

    """

    assert mol is not tmp_mol
    print(mol.inchi)
    print(tmp_mol.inchi)
    assert mol.same(tmp_mol)

    assert mol is not mol2
    assert not mol.same(mol2)


def test_set_position_from_matrix(tmp_mol):

    # The new position matrix just sets all atomic positions to origin.
    new_pos_mat = np.matrix([[0 for x in range(3)] for y in
                            range(tmp_mol.mol.GetNumAtoms())])
    tmp_mol.set_position_from_matrix(new_pos_mat.T)
    for _, atom_coord in tmp_mol.all_atom_coords():
        assert np.allclose(atom_coord, [0, 0, 0], atol=1e-8)


def test_shift(mol):

    s = np.array([10, -20, 5])
    mol2 = mol.shift(s)
    conf = mol2.GetConformer()
    for atom in mol2.GetAtoms():
        atomid = atom.GetIdx()
        pos = conf.GetAtomPosition(atomid)
        should_be = mol.atom_coords(atomid) + s
        assert np.allclose(should_be, pos, atol=1e-8)


def test_update_from_mae(tmp_mol):
    tmp_mol.update_from_mae(join('data',
                                 'molecule',
                                 'molecule.mae'),
                            1)
    assert tmp_mol.max_diameter(0) != tmp_mol.max_diameter(1)
