import itertools as it
import numpy as np
from scipy.spatial.distance import euclidean
import stk
import os
from os.path import join


test_dir = 'molecule_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


def test_all_atom_coords(tmp_amine2):
    """
    Test `all_atom_coords`.

    """

    n_atoms = tmp_amine2.mol.GetNumAtoms()
    new_coords = np.stack([np.arange(n_atoms) for _ in range(3)])
    tmp_amine2.set_position_from_matrix(new_coords)

    for i, (atom, coords) in enumerate(tmp_amine2.all_atom_coords()):
        assert atom == i
        assert all(coords == [i, i, i])

    assert n_atoms == i+1


def test_atom_centroid(amine2):
    assert np.allclose(a=amine2.atom_centroid([0]),
                       b=amine2.atom_coords(0),
                       atol=1e-6)

    atom_ids = [0, 1, 2, 3]
    coords = (amine2.atom_coords(id_) for id_ in atom_ids)
    centroid = sum(coords) / len(atom_ids)
    assert np.allclose(a=amine2.atom_centroid(atom_ids),
                       b=centroid,
                       atol=1e-6)


def test_atom_coords(tmp_amine2):
    """
    Tests `atom_coords`.

    """

    n_atoms = tmp_amine2.mol.GetNumAtoms()
    new_coords = np.stack([np.arange(n_atoms) for _ in range(3)])
    tmp_amine2.set_position_from_matrix(new_coords)

    for i in range(n_atoms):
        coords = tmp_amine2.atom_coords(i)
        assert all(coords == [i, i, i])


def test_atom_distance(amine2):
    """
    Test `atom_distance`.

    """

    # Go through all combinations of atoms in the molecule. Calculate
    # the distance and compare it the distance calculated by the
    # method.
    n_atoms = amine2.mol.GetNumAtoms()
    for atom1_id, atom2_id in it.combinations(range(n_atoms), 2):
        assert (amine2.atom_distance(atom1_id, atom2_id) ==
                euclidean(amine2.atom_coords(atom1_id),
                          amine2.atom_coords(atom2_id)))


def test_atom_symbol(amine2):
    """
    Tests the `atom_symbol` method.

    """

    for atom in amine2.mol.GetAtoms():
        atom_id = atom.GetIdx()
        atom_sym = stk.periodic_table[atom.GetAtomicNum()]
        assert atom_sym == amine2.atom_symbol(atom_id)


def test_center_of_mass(amine2):
    """
    Tests `center_of_mass`.

    """

    # Calculate the center of mass.
    coord_sum = 0
    total_mass = 0
    for atom_id, coord in amine2.all_atom_coords():
        atom_mass = amine2.mol.GetAtomWithIdx(atom_id).GetMass()
        total_mass += atom_mass
        scaled_coord = np.multiply(atom_mass, coord)
        coord_sum = np.add(scaled_coord, coord_sum)

    com = np.divide(coord_sum, total_mass)
    assert np.allclose(amine2.center_of_mass(), com, atol=1e-6)


def test_centroid_functions(tmp_amine2):
    """
    Tests functions related to centroid manipulation of the molecule.

    Functions tested:
        > centroid
        > set_position

    """

    # Save the coordinates of the new centroid.
    new_centroid = tmp_amine2.centroid(0) + np.array([10, 20, 4])
    tmp_amine2.set_position(new_centroid, 0)
    # Check that the centroid is at the desired position.
    assert np.allclose(new_centroid, tmp_amine2.centroid(0), atol=1e-8)


def test_core(amine2):
    for atom in amine2.core().GetAtoms():
        assert atom.GetAtomicNum() != 1
        assert not atom.HasProp('fg')


def test_graph(amine2):
    """
    Tests the output of the `graph` method.

    """

    graph = amine2.graph()
    assert len(graph.nodes()) == amine2.mol.GetNumAtoms()
    assert len(graph.edges()) == amine2.mol.GetNumBonds()


def test_is_core_atom(amine2):
    for atom in amine2.mol.GetAtoms():
        atom_id = atom.GetIdx()
        fg = any(atom_id in fg.atom_ids for fg in amine2.func_groups)
        core = False if fg or atom.GetAtomicNum() == 1 else True
        assert core is amine2.is_core_atom(atom_id)


def test_direction(tmp_amine2):
    direction = tmp_amine2.direction(conformer=1)
    tmp_amine2.set_orientation(direction, [0, 1, 0], conformer=1)
    direction = tmp_amine2.direction(conformer=1)
    assert np.allclose(direction, np.array([0, 1, 0]), atol=1e-4)


def test_position_matrix(amine2):
    position_matrix = amine2.position_matrix()
    for atom_id, atom_coords in amine2.all_atom_coords():
        assert np.array_equal(position_matrix[:, atom_id], atom_coords)


def test_max_diameter(tmp_amine2):
    # Make a position matrix which sets all atoms to the origin except
    # 2 and 13. These should be placed a distance of 100 apart.
    pos_mat = [[0 for x in range(3)] for
               y in range(tmp_amine2.mol.GetNumAtoms())]
    pos_mat[1] = [0, -50, 0]
    pos_mat[3] = [0, 50, 0]
    tmp_amine2.set_position_from_matrix(np.array(pos_mat).T)

    d, id1, id2 = tmp_amine2.max_diameter()
    # Note that it is not exactly 100 because of the Van der Waals
    # radii of the atoms.
    assert d > 100 and d < 105
    assert id1 == 1
    assert id2 == 3


def test_same(amine2,
              tmp_amine2,
              aldehyde2):
    """
    Tests the `same()` method.

    """

    assert amine2 is not tmp_amine2
    assert amine2.same(tmp_amine2)

    assert amine2 is not aldehyde2
    assert not amine2.same(aldehyde2)


def test_set_position_from_matrix(tmp_amine2):
    # The new position matrix just sets all atomic positions to origin.
    new_pos_mat = np.zeros((3, tmp_amine2.mol.GetNumAtoms()))
    tmp_amine2.set_position_from_matrix(new_pos_mat, 0)
    for _, atom_coord in tmp_amine2.all_atom_coords(0):
        assert np.allclose(atom_coord, [0, 0, 0], atol=1e-8)


def test_shift(amine2):

    s = np.array([10, -20, 5])
    mol2 = amine2.shift(s)
    conf = mol2.GetConformer()
    for atom in mol2.GetAtoms():
        atomid = atom.GetIdx()
        pos = conf.GetAtomPosition(atomid)
        should_be = amine2.atom_coords(atomid) + s
        assert np.allclose(should_be, pos, atol=1e-8)


def test_update_from_mae(tmp_amine2, mae_path):
    tmp_amine2.update_from_mae(mae_path, 1)
    assert abs(tmp_amine2.max_diameter(0)[0] -
               tmp_amine2.max_diameter(1)[0]) > 1


def test_update_from_mol(tmp_amine2):
    initial_coords = tmp_amine2.position_matrix()
    new_coords = 4*initial_coords
    # Make a new conformer with new coords.
    tmp_amine2.set_position_from_matrix(new_coords, conformer=1)

    # Write the new conformer to a file.
    path = join(test_dir, 'update_from_mol.mol')
    tmp_amine2.write(
        path=path,
        conformer=1
    )
    # Update the initial conformer with the conformer in the file.
    tmp_amine2.update_from_mol(
        path=path,
        conformer=0
    )

    # The initial conformer should now have the new coords.
    assert np.allclose(
        a=tmp_amine2.position_matrix(conformer=0),
        b=new_coords,
        atol=1e-4
    )
    assert not np.allclose(
        a=tmp_amine2.position_matrix(conformer=0),
        b=initial_coords,
        atol=1e-4
    )


def test_update_from_xyz(tmp_amine2):
    initial_coords = tmp_amine2.position_matrix()
    new_coords = 4*initial_coords
    # Make a new conformer with new coords.
    tmp_amine2.set_position_from_matrix(new_coords, conformer=1)

    # Write the new conformer to a file.
    path = join(test_dir, 'update_from_mol.xyz')
    tmp_amine2.write(
        path=path,
        conformer=1
    )
    # Update the initial conformer with the conformer in the file.
    tmp_amine2.update_from_xyz(
        path=path,
        conformer=0
    )

    # The initial conformer should now have the new coords.
    assert np.allclose(
        a=tmp_amine2.position_matrix(conformer=0),
        b=new_coords,
        atol=1e-4
    )
    assert not np.allclose(
        a=tmp_amine2.position_matrix(conformer=0),
        b=initial_coords,
        atol=1e-4
    )


def test_mol_write(cycle_su):
    atoms = cycle_su.cycle_atoms()

    cycle_su.write(join(test_dir, 'cycle.mol'))
    cycle = stk.StructUnit(join(test_dir, 'cycle.mol'))
    for i in range(cycle.mol.GetNumAtoms()):
        assert np.allclose(
            a=cycle_su.atom_coords(i),
            b=cycle.atom_coords(i),
            atol=1e-3
        )
    for b1 in cycle.mol.GetBonds():
        b2 = cycle_su.mol.GetBondWithIdx(b1.GetIdx())
        assert b1.GetBondType() == b2.GetBondType()
        assert b1.GetBeginAtomIdx() == b2.GetBeginAtomIdx()
        assert b1.GetEndAtomIdx() == b2.GetEndAtomIdx()

    # Write the cycle atoms to a file and create a string for each
    # cycle atom written into valid_lines. Then load the written cycle
    # and for every atom create a string and place it into cycle_lines.
    # Sorted valid and cycle lines should match.
    cycle_su.write(join(test_dir, 'cycle_atoms.mol'), atoms)
    # Note down all atoms written to the file.
    valid_lines = []
    for atom_id in atoms:
        symbol = cycle_su.atom_symbol(atom_id)
        x, y, z = cycle_su.atom_coords(atom_id)
        valid_lines.append(
            f'{symbol} {x:.4f} {y:.4f} {z:.4f}'
        )

    cycle = stk.StructUnit(join(test_dir, 'cycle_atoms.mol'))
    # Note down all atoms read from the file.
    cycle_lines = []
    for atom_id in range(cycle.mol.GetNumAtoms()):
        symbol = cycle.atom_symbol(atom_id)
        x, y, z = cycle.atom_coords(atom_id)
        cycle_lines.append(
            f'{symbol} {x:.4f} {y:.4f} {z:.4f}'
        )

    # Written and read atoms should match.
    assert sorted(valid_lines) == sorted(cycle_lines)


def test_xyz_write(cycle_su):
    atoms = cycle_su.cycle_atoms()
    cycle_su.write(join(test_dir, 'cycle.xyz'))

    with open(join(test_dir, 'cycle.xyz'), 'r') as f:
        xyz_content = f.read()

    # Check that every atom in the molecule has a line in the file.
    xyz_lines = set(xyz_content.split('\n'))
    for atom_id in range(cycle_su.mol.GetNumAtoms()):
        x, y, z = cycle_su.atom_coords(atom_id)
        symbol = cycle_su.atom_symbol(atom_id)
        assert f'{symbol} {x:f} {y:f} {z:f}' in xyz_lines

    assert len(xyz_content.split('\n'))-3 == cycle_su.mol.GetNumAtoms()

    # Test writing of specific atoms.
    cycle_su.write(join(test_dir, 'cycle_atoms.xyz'), atoms)
    with open(join(test_dir, 'cycle_atoms.xyz'), 'r') as f:
        xyz_content = f.read()

    # Check that every cycle atom has a line in the file.
    xyz_lines = set(xyz_content.split('\n'))
    for atom_id in atoms:
        x, y, z = cycle_su.atom_coords(atom_id)
        symbol = cycle_su.atom_symbol(atom_id)
        assert f'{symbol} {x:f} {y:f} {z:f}' in xyz_lines

    # Check that only cycle atoms were written.
    assert len(xyz_content.split('\n'))-3 == len(atoms)


def test_pdb_write(cycle_su):
    cycle_su.write(join(test_dir, 'cycle.pdb'))
    cycle = stk.StructUnit(join(test_dir, 'cycle.pdb'))

    # Make sure the position matrices are basically the same.
    assert np.allclose(
        a=cycle_su.position_matrix(),
        b=cycle.position_matrix(),
        atol=1e-5
    )

    # Make sure the connectivity is the same.
    bonds1 = sorted(
        sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
        for b in cycle_su.mol.GetBonds()
    )

    bonds2 = sorted(
        sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
        for b in cycle.mol.GetBonds()
    )
    assert bonds1 == bonds2

    # Test writing of specific atoms only.
    atoms = cycle_su.cycle_atoms()
    cycle_su.write(join(test_dir, 'cycle_atoms.pdb'), atoms)
    cycle = stk.StructUnit(join(test_dir, 'cycle_atoms.pdb'))

    # Make sure the position matrices are basically the same.
    assert np.allclose(
        a=cycle_su.position_matrix()[:, atoms],
        b=cycle.position_matrix(),
        atol=1e-5
    )

    # Make sure the bonds which need to exist exist.
    # Atom ids will change, so compare positions.
    atoms = set(atoms)
    bonds1 = []
    for bond in cycle_su.mol.GetBonds():
        a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a1 in atoms and a2 in atoms:
            bonds1.append(sorted([
                *cycle_su.atom_coords(a1),
                *cycle_su.atom_coords(a2)
            ]))
    bonds1 = np.array(bonds1)

    bonds2 = []
    for bond in cycle.mol.GetBonds():
        a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bonds2.append(sorted([
            *cycle.atom_coords(a1),
            *cycle.atom_coords(a2)
        ]))
    bonds2 = np.array(bonds2)

    assert np.allclose(
        a=bonds1,
        b=bonds2,
        atol=1e-5
    )
