import os
import pytest
import numpy as np
import itertools as it
import stk


def test_apply_displacement(molecule, displacement):
    before = molecule.get_position_matrix()
    molecule.apply_displacement(displacement)
    assert np.allclose(
        a=before+displacement,
        b=molecule.get_position_matrix(),
        atol=1e-32,
    )


class TestApplyRotationAboutAxis:
    def _rotational_space_positions(self, molecule, axis, origin):
        """
        Get the atomic coordinates on the plane of the rotation.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            The molecule being rotated.

        axis : :class:`numpy.ndarray`
            The axis about which the rotation happens.

        origin : :class:`numpy.ndarray`
            The origin about which the rotation happens.

        Returns
        -------
        :class:`numpy.ndarray`
            An ``[n, 3]`` of atomic positions of `molecule`, projected
            onto the plane about which the rotation happens. The
            `axis` is the normal to this plane.

        """

        axis_matrix = np.repeat([axis], len(molecule.atoms), 0).T
        positions = molecule.get_position_matrix() - origin
        return positions - (axis_matrix * (positions @ axis)).T

    def test(
        self,
        molecule,
        angle,
        axis,
        origin
    ):
        before = self._rotational_space_positions(
            molecule=molecule,
            axis=axis,
            origin=origin,
        )
        molecule.apply_rotation_about_axis(angle, axis, origin)
        after = self._rotational_space_positions(
            molecule=molecule,
            axis=axis,
            origin=origin,
        )

        for atom_id in range(len(molecule.atoms)):
            applied_rotation = stk.vector_angle(
                vector1=before[atom_id],
                vector2=after[atom_id],
            )
            assert abs(abs(angle) - applied_rotation) < 1e-13


def test_get_atom_positions(molecule, get_atom_ids_no_fail):
    position_matrix = molecule.get_position_matrix()
    atom_ids = get_atom_ids_no_fail(molecule)
    positions = it.zip_longest(
        range(len(molecule.atoms)) if atom_ids is None else atom_ids,
        molecule.get_atom_positions(get_atom_ids_no_fail(molecule)),
    )
    for atom_id, position in positions:
        assert np.all(np.equal(position, position_matrix[atom_id]))


def test_get_atom_distance(molecule):
    position_matrix = molecule.get_position_matrix()
    positions_1 = np.repeat([position_matrix], len(molecule.atoms), 0)
    positions_2 = positions_1.swapaxes(0, 1)
    distance_matrix = np.linalg.norm(positions_1 - positions_2, axis=2)

    atom_ids = range(len(molecule.atoms))
    for atom1, atom2 in it.product(atom_ids, atom_ids):
        true_distance = distance_matrix[atom1, atom2]
        distance = molecule.get_atom_distance(atom1, atom2)
        assert abs(true_distance - distance) < 1e-13


def test_get_cached_mol():
    bb1 = stk.BuildingBlock('NCCN')
    with pytest.raises(KeyError):
        stk.BuildingBlock.get_cached_mol(bb1.get_identity_key())

    bb2 = stk.BuildingBlock('NCCN', use_cache=True)
    assert (
        stk.BuildingBlock.get_cached_mol(bb2.get_identity_key()) is bb2
    )


def test_get_center_of_mass(molecule, get_atom_ids):
    atom_ids = get_atom_ids(molecule)
    center_of_mass = molecule.get_center_of_mass(atom_ids),

    if atom_ids is None:
        atom_ids = range(len(molecule.atoms))
    else:
        atom_ids = list(get_atom_ids(molecule))

    valid_atom_ids = set(atom_ids)
    atoms = filter(
        lambda atom: atom.id in valid_atom_ids,
        molecule.atoms
    )
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
        atom_ids = range(len(molecule.atoms))
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
        bb.set_position_matrix(
            np.array([[i, 0, 0] for i in range(len(bb.atoms))])
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
        bb.set_position_matrix(coords)

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
        num_atoms = len(molecule.atoms)
        coords = np.array([[i, 0, 0] for i in range(num_atoms)])
        molecule.set_position_matrix(coords)
        return molecule, None, num_atoms-1

    def case2(atom_ids, maximum_diameter):
        molecule = stk.BuildingBlock('NCCN')
        num_atoms = len(molecule.atoms)
        coords = np.zeros((num_atoms, 3))
        coords[[1]] = [0, -50, 0]
        coords[[9]] = [0, 50, 0]
        molecule.set_position_matrix(coords)
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
        molecule.set_position_matrix(coords)
        return molecule, atom_ids, normal

    def case2(atom_ids, normal):
        molecule = stk.BuildingBlock('NCCN')
        coords = molecule.get_position_matrix()
        coords[:, 2] = 0
        molecule.set_position_matrix(coords)
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
        return molecule, np.zeros((len(molecule.atoms), 3))

    def case2():
        molecule = stk.BuildingBlock('NCCN')
        return molecule, np.ones((len(molecule.atoms), 3))

    @pytest.mark.parametrize(
        'molecule,position_matrix',
        [
            case1(),
            case2(),
        ],
    )
    def test(self, molecule, position_matrix):
        molecule.set_position_matrix(position_matrix)
        assert np.all(np.equal(
            molecule.get_position_matrix(),
            position_matrix,
        ))


class TestUpdateFromFile1:
    @pytest.fixture(params=[
        'molecule.mol',
        'molecule.xyz',
    ])
    def path(self, request, tmpdir):
        return os.path.join(tmpdir, request.param)

    def case1():
        conformer1 = stk.BuildingBlock('NCCN')
        conformer2 = stk.BuildingBlock('NCCN')
        conformer2.set_position_matrix(
            position_matrix=np.zeros((len(conformer2.atoms), 3)),
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
            position_matrix=np.zeros((len(conformer2.atoms), 3)),
        )
        return conformer1, conformer2

    @pytest.mark.parametrize(
        'get_conformers',
        [
            case1,
            case2,
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
        conformer1.update_from_file(path)

        assert np.allclose(
            a=conformer1.get_position_matrix(),
            b=conformer2.get_position_matrix(),
            atol=1e-4,
        )


def test_update_from_rdkit_mol(molecule):
    before = molecule.get_position_matrix()
    rdkit_molecule = molecule.to_rdkit_mol()
    conformer = rdkit_molecule.GetConformer()
    for atom_id, position in enumerate(conformer.GetPositions()):
        conformer.SetAtomPosition(atom_id, 0.5*position)

    molecule.update_from_rdkit_mol(rdkit_molecule)
    after = molecule.get_position_matrix()
    assert np.allclose(conformer.GetPositions(), after, 1e-32)
    assert not np.allclose(before, after, 1e-32)


def test_to_rdkit_mol(molecule):
    rdkit_molecule = molecule.to_rdkit_mol()
    assert rdkit_molecule.GetNumConformers() == 1

    assert rdkit_molecule.GetNumAtoms() == len(molecule.atoms)
    atoms = zip(molecule.atoms, rdkit_molecule.GetAtoms())
    for atom, rdkit_atom in atoms:
        assert atom.charge == rdkit_atom.GetFormalCharge()
        assert atom.atomic_number == rdkit_atom.GetAtomicNum()
        assert atom.mass == rdkit_atom.GetMass()

    assert len(molecule.bonds) == rdkit_molecule.GetNumBonds()
    bonds = zip(molecule.bonds, rdkit_molecule.GetBonds())
    for bond, rdkit_bond in bonds:
        assert bond.order == rdkit_bond.GetBondTypeAsDouble()
        assert bond.atom1.id == rdkit_bond.GetBeginAtomIdx()
        assert bond.atom2.id == rdkit_bond.GetEndAtomIdx()
