import os
import pytest
import numpy as np
import itertools as it
import stk


def test_get_atomic_positions(molecule, get_atom_ids_no_fail):
    position_matrix = molecule.get_position_matrix()
    atom_ids = get_atom_ids_no_fail(molecule)
    if atom_ids is None:
        atom_ids = range(molecule.get_num_atoms())

    positions = it.zip_longest(
        atom_ids,
        molecule.get_atomic_positions(get_atom_ids_no_fail(molecule)),
    )
    for atom_id, position in positions:
        assert np.all(np.equal(position, position_matrix[atom_id]))


def test_get_atomic_distance(molecule):
    position_matrix = molecule.get_position_matrix()
    positions_1 = np.repeat(
        a=[position_matrix],
        repeats=molecule.get_num_atoms(),
        axis=0,
    )
    positions_2 = positions_1.swapaxes(0, 1)
    distance_matrix = np.linalg.norm(positions_1 - positions_2, axis=2)

    atom_ids = range(molecule.get_num_atoms())
    for atom1, atom2 in it.product(atom_ids, atom_ids):
        true_distance = distance_matrix[atom1, atom2]
        distance = molecule.get_atomic_distance(atom1, atom2)
        assert abs(true_distance - distance) < 1e-13


def test_get_center_of_mass(molecule, get_atom_ids):
    atom_ids = get_atom_ids(molecule)
    center_of_mass = molecule.get_center_of_mass(atom_ids),

    if atom_ids is None:
        atom_ids = range(molecule.get_num_atoms())
    else:
        atom_ids = list(get_atom_ids(molecule))

    atoms = molecule.get_atoms(atom_ids=atom_ids)
    masses = [[atom.mass] for atom in atoms]
    true_center_of_mass = np.divide(
        np.sum(
            a=masses*molecule.get_position_matrix()[atom_ids, :],
            axis=0,
        ),
        sum(mass for mass, in masses),
    )
    assert np.allclose(
        a=true_center_of_mass,
        b=center_of_mass,
        atol=1e-32,
    )


def test_get_centroid(molecule, get_atom_ids):
    atom_ids = get_atom_ids(molecule)
    centroid = molecule.get_centroid(atom_ids)

    if atom_ids is None:
        atom_ids = range(molecule.get_num_atoms())
    else:
        atom_ids = list(get_atom_ids(molecule))

    true_centroid = np.divide(
        np.sum(
            a=molecule.get_position_matrix()[atom_ids, :],
            axis=0
        ),
        len(atom_ids),
    )
    assert np.allclose(
        a=true_centroid,
        b=centroid,
        atol=1e-32,
    )


class TestGetDirection1:
    def case1():
        bb = stk.BuildingBlock('NCCN')
        bb = bb.with_position_matrix(
            np.array([[i, 0, 0] for i in range(bb.get_num_atoms())])
        )
        return bb, [1, 0, 0]

    @pytest.mark.parametrize(
        'molecule,direction',
        [
            case1(),
        ],
    )
    def test(self, molecule, get_atom_ids, direction):
        atom_ids = get_atom_ids(molecule)
        assert np.allclose(
            a=molecule.get_direction(atom_ids),
            b=direction,
            atol=1e-32,
        )


class TestGetDirection2:
    def case1():
        bb = stk.BuildingBlock('NCCN')
        atom_ids = [1, 3]

        coords = bb.get_position_matrix()
        coords[atom_ids] = [[1, 1, 1], [3, 3, 3]]
        bb = bb.with_position_matrix(coords)

        return bb, atom_ids, [1/np.sqrt(3)]*3

    @pytest.mark.parametrize(
        'molecule,atom_ids,direction',
        [
            case1(),
        ],
    )
    def test(self, molecule, atom_ids, direction):
        assert np.allclose(
            a=molecule.get_direction(atom_ids),
            b=direction,
            atol=1e-32,
        )


class TestGetMaximumDiameter:
    def case1():
        molecule = stk.BuildingBlock('NCCN')
        coords = np.array([
            [i, 0, 0] for i in range(molecule.get_num_atoms())
        ])
        molecule = molecule.with_position_matrix(coords)
        return molecule, None, molecule.get_num_atoms()-1

    def case2(atom_ids, maximum_diameter):
        molecule = stk.BuildingBlock('NCCN')
        coords = np.zeros((molecule.get_num_atoms(), 3))
        coords[[1]] = [0, -50, 0]
        coords[[9]] = [0, 50, 0]
        molecule = molecule.with_position_matrix(coords)
        return molecule, atom_ids, maximum_diameter

    @pytest.mark.parametrize(
        'molecule,atom_ids,maximum_diameter',
        [
            case1(),
            case2(atom_ids=None, maximum_diameter=100),
            case2(atom_ids=(1, 9), maximum_diameter=100),
            case2(atom_ids=(1, 0), maximum_diameter=50),
            case2(atom_ids=(0, 9), maximum_diameter=50),
            case2(atom_ids=(0, 2, 3, 4), maximum_diameter=0),
        ],
    )
    def test(
        self,
        molecule,
        atom_ids,
        maximum_diameter
    ):
        assert (
            molecule.get_maximum_diameter(atom_ids) == maximum_diameter
        )


class TestGetPlaneNormal:
    def case1(atom_ids, normal):
        molecule = stk.BuildingBlock('NCCN')
        coords = molecule.get_position_matrix()
        coords[[1, 9], 2] = 0
        molecule = molecule.with_position_matrix(coords)
        return molecule, atom_ids, normal

    def case2(atom_ids, normal):
        molecule = stk.BuildingBlock('NCCN')
        coords = molecule.get_position_matrix()
        coords[:, 2] = 0
        molecule = molecule.with_position_matrix(coords)
        return molecule, atom_ids, normal

    @pytest.mark.parametrize(
        'molecule,atom_ids,normal',
        [
            case1(atom_ids=(1, 9), normal=[0, 0, 1]),
            case2(atom_ids=None, normal=[0, 0, 1]),
            case2(atom_ids=(1, 9), normal=[0, 0, 1]),
        ],
    )
    def test_1(self, molecule, atom_ids, normal):
        assert np.all(np.equal(
            molecule.get_plane_normal(atom_ids),
            normal,
        ))

    @pytest.mark.parametrize(
        'molecule,atom_ids,normal',
        [
            case1(atom_ids=None, normal=[0, 0, 1])
        ],
    )
    def test_2(self, molecule, atom_ids, normal):
        assert not np.all(np.equal(
            molecule.get_plane_normal(atom_ids),
            normal,
        ))


class TestPositionMatrix:
    def case1():
        molecule = stk.BuildingBlock('NCCN')
        return molecule, np.zeros((molecule.get_num_atoms(), 3))

    def case2():
        molecule = stk.BuildingBlock('NCCN')
        return molecule, np.ones((molecule.get_num_atoms(), 3))

    @pytest.mark.parametrize(
        'molecule,position_matrix',
        [
            case1(),
            case2(),
        ],
    )
    def test(self, molecule, position_matrix):
        # Keep clone to test for immutability.
        clone = molecule.clone()

        new = molecule.with_position_matrix(position_matrix)
        assert np.all(np.equal(
            new.get_position_matrix(),
            position_matrix,
        ))
        _test_unchanged(molecule, clone)


def test_with_centroid(molecule, get_atom_ids):
    # Keep clone to test for immutability.
    clone = molecule.clone()

    new = molecule.with_centroid(
        position=[1, 2, 3],
        atom_ids=get_atom_ids(molecule),
    )
    assert np.allclose(
        a=new.get_centroid(atom_ids=get_atom_ids(molecule)),
        b=[1, 2, 3],
        atol=1e-32,
    )
    _test_unchanged(clone, molecule)


class TestWithStructureFromFile1:
    @pytest.fixture(params=[
        'molecule.mol',
        'molecule.xyz',
    ])
    def path(self, request, tmpdir):
        return os.path.join(tmpdir, request.param)

    def case1():
        conformer1 = stk.BuildingBlock('NCCN')
        conformer2 = stk.BuildingBlock('NCCN')
        conformer2 = conformer2.with_position_matrix(
            position_matrix=np.zeros((conformer2.get_num_atoms(), 3)),
        )
        return conformer1, conformer2

    def case2():
        bb = stk.BuildingBlock('BrCCBr', ['bromine'])
        topology_graph = stk.polymer.Linear('A', 3)
        conformer1 = stk.ConstructedMolecule(
            building_blocks=[bb],
            topology_graph=topology_graph,
        )
        conformer2 = stk.ConstructedMolecule(
            building_blocks=[bb],
            topology_graph=topology_graph,
        )
        conformer2.set_position_matrix(
            position_matrix=np.zeros((conformer2.get_num_atoms(), 3)),
        )
        return conformer1, conformer2

    @pytest.mark.parametrize(
        'get_conformers',
        [
            case1,
            # case2,
        ],
    )
    def test(self, get_conformers, path):
        conformer1, conformer2 = get_conformers()
        assert not np.allclose(
            a=conformer1.get_position_matrix(),
            b=conformer2.get_position_matrix(),
            atol=1e-4,
        )

        conformer2.write(path)
        # Keep a clone for immutability testing.
        clone = conformer1.clone()
        new = conformer1.with_structure_from_file(path)

        assert np.allclose(
            a=new.get_position_matrix(),
            b=conformer2.get_position_matrix(),
            atol=1e-4,
        )
        _test_unchanged(clone, conformer1)


class TestWithStructureFromFile2:
    @pytest.fixture
    def path(self, datadir, request):
        return str(datadir / request.param)

    def case1():
        return stk.BuildingBlock('NCCN'), 'NCCN.mae'

    @pytest.mark.parametrize(
        'molecule,path',
        [
            case1(),
        ],
        indirect=['path'],
    )
    def test(self, molecule, path):
        # Keep clone for immutability testing.
        clone = molecule.clone()
        new = molecule.with_structure_from_file(path)
        size_diff = abs(
            molecule.get_maximum_diameter()
            - new.get_maximum_diameter()
        )
        assert size_diff > 1
        _test_unchanged(clone, molecule)


def test_to_rdkit_mol(molecule):
    rdkit_molecule = molecule.to_rdkit_mol()
    assert rdkit_molecule.GetNumConformers() == 1

    assert rdkit_molecule.GetNumAtoms() == molecule.get_num_atoms()
    atoms = zip(molecule.get_atoms(), rdkit_molecule.GetAtoms())
    for atom, rdkit_atom in atoms:
        assert atom.charge == rdkit_atom.GetFormalCharge()
        assert atom.atomic_number == rdkit_atom.GetAtomicNum()
        assert atom.mass == rdkit_atom.GetMass()

    assert molecule.get_num_bonds() == rdkit_molecule.GetNumBonds()
    bonds = zip(molecule.get_bonds(), rdkit_molecule.GetBonds())
    for bond, rdkit_bond in bonds:
        assert bond.order == rdkit_bond.GetBondTypeAsDouble()
        assert bond.atom1.id == rdkit_bond.GetBeginAtomIdx()
        assert bond.atom2.id == rdkit_bond.GetEndAtomIdx()


def is_equivalent_atom(atom1, atom2):
    return (
        atom1 is not atom2
        and atom1.id == atom2.id
        and atom1.charge == atom2.charge
        and atom1.__class__ is atom2.__class__
    )


def is_equivalent_bond(bond1, bond2):
    return (
        bond1 is not bond2
        and bond1.__class__ is bond2.__class__
        and bond1.order == bond2.order
        and bond1.periodicity == bond2.periodicity
        and is_equivalent_atom(bond1.atom1, bond2.atom1)
        and is_equivalent_atom(bond1.atom2, bond2.atom2)
    )


def test_clone(molecule):
    clone = molecule.clone()
    atoms = it.zip_longest(clone.get_atoms(), molecule.get_atoms())
    for a1, a2 in atoms:
        assert is_equivalent_atom(a1, a2)

    bonds = it.zip_longest(clone.get_bonds(), molecule.get_bonds())
    for b1, b2 in bonds:
        assert is_equivalent_bond(b1, b2)
