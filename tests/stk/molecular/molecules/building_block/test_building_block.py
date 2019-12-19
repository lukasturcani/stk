import numpy as np
import itertools as it
import stk
import pytest
import rdkit.Chem.AllChem as rdkit


class _FunctionalGroup:
    def __init__(self, atom_ids, bonder_ids, deleter_ids):
        self.atom_ids = atom_ids
        self.bonder_ids = bonder_ids
        self.deleter_ids = deleter_ids


class TestInitFromRdkitMol:
    def case1():
        rdkit_molecule = rdkit.AddHs(rdkit.MolFromSmiles('NCCN'))
        rdkit.EmbedMolecule(rdkit_molecule, rdkit.ETKDGv2())

        functional_groups = ['amine']

        expected_functional_groups = [
            _FunctionalGroup(
                atom_ids=(0, 4, 5),
                bonder_ids=(0, ),
                deleter_ids=(4, 5),
            ),
            _FunctionalGroup(
                atom_ids=(3, 10, 11),
                bonder_ids=(3, ),
                deleter_ids=(10, 11),
            ),
        ]
        return (
            rdkit_molecule,
            functional_groups,
            expected_functional_groups,
        )

    def case2():
        rdkit_molecule = rdkit.AddHs(rdkit.MolFromSmiles('NCCN'))
        rdkit.EmbedMolecule(rdkit_molecule, rdkit.ETKDGv2())
        return rdkit_molecule, None, []

    def case3():
        rdkit_molecule = rdkit.AddHs(rdkit.MolFromSmiles('NCCN'))
        rdkit.EmbedMolecule(rdkit_molecule, rdkit.ETKDGv2())
        return rdkit_molecule, ['aldehyde'], []

    @pytest.mark.parametrize(
        argnames=(
            'rdkit_molecule',
            'functional_groups',
            'expected_functional_groups',
        ),
        argvalues=[
            case1(),
            case2(),
            case3(),
        ],
    )
    def test(
        self,
        rdkit_molecule,
        functional_groups,
        expected_functional_groups,
    ):
        molecule = stk.BuildingBlock.init_from_rdkit_mol(
            mol=rdkit_molecule,
            functional_groups=functional_groups,
        )
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

        num_func_groups = len(molecule.func_groups)
        assert num_func_groups == len(expected_functional_groups)
        func_groups = zip(
            molecule.func_groups,
            expected_functional_groups,
        )
        for fg, expected_fg in func_groups:
            assert tuple(fg.get_atom_ids()) == expected_fg.atom_ids
            assert tuple(fg.get_bonder_ids()) == expected_fg.bonder_ids
            assert (
                tuple(fg.get_deleter_ids()) == expected_fg.deleter_ids
            )


def is_equivalent_atom(atom1, atom2):
    return (
        atom1.id == atom2.id
        and atom1.charge == atom2.charge
        and atom1.__class__ is atom2.__class__
    )


def is_equivalent_bond(bond1, bond2):
    return (
        bond1.__class__ is bond2.__class__
        and bond1.order == bond2.order
        and is_equivalent_atom(bond1.atom1, bond2.atom1)
        and is_equivalent_atom(bond1.atom2, bond2.atom2)
        and bond1.periodicity == bond2.periodicity
    )


def is_equivalent_fg(fg1, fg2):
    equivalent_atoms = all(
        id1 == id2
        for id1, id2
        in it.zip_longest(fg1.get_atom_ids(), fg2.get_atom_ids())
    )
    equivalent_bonders = all(
        id1 == id2
        for id1, id2
        in it.zip_longest(fg1.get_bonder_ids(), fg2.get_bonder_ids())
    )
    equivalent_deleters = all(
        id1 == id2
        for id1, id2
        in it.zip_longest(fg1.get_deleter_ids(), fg2.get_deleter_ids())
    )
    return (
        fg1.fg_type is fg2.fg_type
        and equivalent_atoms
        and equivalent_bonders
        and equivalent_deleters
    )


def is_equivalent_building_block(building_block1, building_block2):
    atoms = it.zip_longest(
        building_block1.atoms,
        building_block2.atoms,
    )
    equivalent_atoms = all(
        is_equivalent_atom(a1, a2) for a1, a2 in atoms
    )

    bonds = it.zip_longest(
        building_block1.bonds,
        building_block2.bonds,
    )
    equivalent_bonds = all(
        is_equivalent_bond(b1, b2) for b1, b2 in bonds
    )

    fgs = it.zip_longest(
        building_block1.func_groups,
        building_block2.func_groups,
    )
    equivalent_fgs = all(
        is_equivalent_fg(fg1, fg2) for fg1, fg2 in fgs
    )
    return equivalent_atoms and equivalent_bonds and equivalent_fgs


class TestInitFromFile:
    @pytest.fixture(
        params=[
            'building_block.mol',
            'building_block.pdb',
        ],
    )
    def filename(self, request):
        return request.param

    def test(self, tmpdir, filename, building_block):
        path = str(tmpdir / filename)
        building_block.write(path)

        loaded = stk.BuildingBlock.init_from_file(
            path=path,
            functional_groups={
                fg.fg_type.name for fg in building_block.func_groups
            },
        )

        atoms = it.zip_longest(building_block.atoms, loaded.atoms)
        for a1, a2 in atoms:
            assert is_equivalent_atom(a1, a2)

        bonds = it.zip_longest(building_block.bonds, loaded.bonds)
        for b1, b2 in bonds:
            assert is_equivalent_bond(b1, b2)

        fgs = it.zip_longest(
            building_block.func_groups,
            loaded.func_groups
        )
        for fg1, fg2 in fgs:
            assert is_equivalent_fg(fg1, fg2)


def test_clone(building_block):
    clone = building_block.clone()
    fgs = it.zip_longest(clone.func_groups, building_block.func_groups)
    for fg1, fg2 in fgs:
        assert is_equivalent_fg(fg1, fg2)


class TestInitFromMolecule:
    def case1():
        functional_groups = ['amine']
        building_block = stk.BuildingBlock(
            smiles='NCCN',
            functional_groups=functional_groups,
        )
        return (
            building_block,
            functional_groups,
            building_block.atoms,
            building_block.bonds,
            building_block.func_groups,
        )

    def case2():
        functional_groups = ['bromine']
        molecule = stk.ConstructedMolecule(
            building_blocks=[
                stk.BuildingBlock('BrCCBr', functional_groups),
            ],
            topology_graph=stk.polymer.Linear('AA', 1),
        )
        expected_functional_groups = (
            stk.FunctionalGroup(
                atoms=(stk.C(1), stk.Br(2)),
                bonders=(stk.C(1), ),
                deleters=(stk.Br(2), ),
                fg_type=stk.functional_groups.fg_types['bromine'],
            ),
            stk.FunctionalGroup(
                atoms=(stk.C(8), stk.Br(7)),
                bonders=(stk.C(8), ),
                deleters=(stk.Br(7), ),
                fg_type=stk.functional_groups.fg_types['bromine'],
            ),
        )
        return (
            molecule,
            functional_groups,
            molecule.atoms,
            molecule.bonds,
            expected_functional_groups,
        )

    @pytest.mark.parametrize(
        argnames=(
            'molecule',
            'functional_groups',
            'expected_atoms',
            'expected_bonds',
            'expected_functional_groups',
        ),
        argvalues=(
            case1(),
            case2(),
        ),
    )
    def test(
        self,
        molecule,
        functional_groups,
        expected_atoms,
        expected_bonds,
        expected_functional_groups,
    ):
        building_block = stk.BuildingBlock.init_from_molecule(
            mol=molecule,
            functional_groups=functional_groups,
        )
        atoms = it.zip_longest(
            building_block.atoms,
            expected_atoms,
        )
        for a1, a2 in atoms:
            assert is_equivalent_atom(a1, a2)

        bonds = it.zip_longest(
            building_block.bonds,
            expected_bonds,
        )
        for b1, b2 in bonds:
            assert is_equivalent_bond(b1, b2)

        fgs = it.zip_longest(
            building_block.func_groups,
            expected_functional_groups,
        )
        for fg1, fg2 in fgs:
            assert is_equivalent_fg(fg1, fg2)


class TestInitFromRandomFile:
    @pytest.fixture
    def file_glob(self, datadir, request):
        return str(datadir / request.param)

    def case1():
        functional_groups = ['amine']
        building_block = stk.BuildingBlock(
            smiles='NCCCN',
            functional_groups=functional_groups,
        )
        return 'neutral.mol', functional_groups, building_block

    def case2():
        functional_groups = ['amine']
        building_block = stk.BuildingBlock(
            smiles='NC[C-]CN',
            functional_groups=functional_groups,
        )
        return 'negative_carbon.mol', functional_groups, building_block

    def case3():
        functional_groups = ['amine']
        building_block = stk.BuildingBlock(
            smiles='[N-]CCCN',
            functional_groups=functional_groups,
        )
        return (
            'negative_nitrogen.mol', functional_groups, building_block
        )

    @pytest.mark.parametrize(
        argnames=(
            'file_glob',
            'functional_groups',
            'expected_building_block',
        ),
        argvalues=(
            case1(),
            case2(),
            case3(),
        ),
        indirect=['file_glob'],
    )
    def test(
        self,
        file_glob,
        functional_groups,
        expected_building_block,
    ):
        building_block = stk.BuildingBlock.init_from_random_file(
            file_glob=file_glob,
            functional_groups=functional_groups,
        )
        assert is_equivalent_building_block(
            building_block1=building_block,
            building_block2=expected_building_block,
        )


class TestInitFromSmiles:
    def case1():
        expected_atoms = (
            stk.N(0),
            stk.C(1),
            stk.C(2),
            stk.N(3),
            stk.H(4),
            stk.H(5),
            stk.H(6),
            stk.H(7),
            stk.H(8),
            stk.H(9),
            stk.H(10),
            stk.H(11),
        )
        expected_bonds = (
            stk.Bond(stk.N(0), stk.C(1), 1),
            stk.Bond(stk.C(1), stk.C(2), 1),
            stk.Bond(stk.C(2), stk.N(3), 1),
            stk.Bond(stk.N(0), stk.H(4), 1),
            stk.Bond(stk.N(0), stk.H(5), 1),
            stk.Bond(stk.C(1), stk.H(6), 1),
            stk.Bond(stk.C(1), stk.H(7), 1),
            stk.Bond(stk.C(2), stk.H(8), 1),
            stk.Bond(stk.C(2), stk.H(9), 1),
            stk.Bond(stk.N(3), stk.H(10), 1),
            stk.Bond(stk.N(3), stk.H(11), 1),
        )
        expected_functional_groups = (
            stk.FunctionalGroup(
                atoms=(stk.N(0), stk.H(4), stk.H(5)),
                bonders=(stk.N(0),),
                deleters=(stk.H(4), stk.H(5)),
                fg_type=stk.functional_groups.fg_types['amine'],
            ),
            stk.FunctionalGroup(
                atoms=(stk.N(3), stk.H(10), stk.H(11)),
                bonders=(stk.N(3), ),
                deleters=(stk.H(10), stk.H(11)),
                fg_type=stk.functional_groups.fg_types['amine'],
            ),
        )
        return (
            'NCCN',
            ['amine'],
            expected_atoms,
            expected_bonds,
            expected_functional_groups,
        )

    @pytest.mark.parametrize(
        argnames=(
            'smiles',
            'functional_groups',
            'expected_atoms',
            'expected_bonds',
            'expected_functional_groups',
        ),
        argvalues=(
            case1(),
        ),
    )
    def test(
        self,
        smiles,
        functional_groups,
        expected_atoms,
        expected_bonds,
        expected_functional_groups,
    ):
        building_block = stk.BuildingBlock(smiles, functional_groups)
        atoms = it.zip_longest(
            building_block.atoms,
            expected_atoms,
        )
        for a1, a2 in atoms:
            assert is_equivalent_atom(a1, a2)

        bonds = it.zip_longest(
            building_block.bonds,
            expected_bonds,
        )
        for b1, b2 in bonds:
            assert is_equivalent_bond(b1, b2)

        fgs = it.zip_longest(
            building_block.func_groups,
            expected_functional_groups,
        )
        for fg1, fg2 in fgs:
            assert is_equivalent_fg(fg1, fg2)


def test_get_bonder_ids(building_block, get_fg_ids):
    fg_ids = get_fg_ids(building_block)
    if fg_ids is None:
        fg_ids = range(len(building_block.func_groups))
    fg_ids = list(fg_ids)
    expected_bonder_ids = (
        bid
        for fg_id in fg_ids
        for bid in building_block.func_groups[fg_id].get_bonder_ids()
    )
    ids = building_block.get_bonder_ids(
        fg_ids=get_fg_ids(building_block),
    )
    bonder_ids = it.zip_longest(ids, expected_bonder_ids)
    for bonder_id, expected_bonder_id in bonder_ids:
        assert bonder_id == expected_bonder_id


class TestGetBonderCentroids:
    def case1():
        building_block = stk.BuildingBlock('BrCCBr', ['bromine'])
        position_matrix = np.zeros((len(building_block.atoms), 3))
        bonder1 = [1.0, 0., 0.]
        bonder2 = [10., 0., 0.]
        position_matrix[1, :] = bonder1
        position_matrix[2, :] = bonder2
        building_block.set_position_matrix(position_matrix)
        fg_ids = (0, 1)
        return building_block, fg_ids, [bonder1, bonder2]

    @pytest.mark.parametrize(
        argnames=(
            'building_block',
            'fg_ids',
            'expected_bonder_centroids',
        ),
        argvalues=(
            case1(),
        )
    )
    def test(
        self,
        building_block,
        fg_ids,
        expected_bonder_centroids,
    ):
        bonder_centroids = building_block.get_bonder_centroids(
            fg_ids=fg_ids,
        )
        bonder_centroids = it.zip_longest(
            bonder_centroids,
            expected_bonder_centroids,
        )
        for centroid, expected_bonder_centroid in bonder_centroids:
            assert np.allclose(
                a=centroid,
                b=expected_bonder_centroid,
                atol=1e-32,
            )


class _BonderPlanePoint:
    def __init__(self, x, y, z, distance):
        self.x = x
        self.y = y
        self.z = z
        self.distance = distance


class TestGetBonderPlane:

    def case1():
        building_block = stk.BuildingBlock(
            smiles='BrCC(CBr)(CBr)CBr',
            functional_groups=['bromine'],
        )
        bonder_ids = list(building_block.get_bonder_ids())
        coords = building_block.get_position_matrix()
        coords[bonder_ids[0]] = [1, 1, 0]
        coords[bonder_ids[1]] = [0, 0, 0.5]
        coords[bonder_ids[2]] = [0, 0, -0.5]
        coords[bonder_ids[3]] = [1, -1, 0]
        building_block.set_position_matrix(coords)

        points = (
            _BonderPlanePoint(1, 1, 0, 0),
            _BonderPlanePoint(0, 0, 0.5, 0.5),
            _BonderPlanePoint(0, 0, -0.5, 0.5),
            _BonderPlanePoint(1, -1, 0, 0),
        )

        return building_block, None, points

    @pytest.mark.parametrize(
        argnames=(
            'building_block',
            'fg_ids',
            'points',
        ),
        argvalues=(
            case1(),
        ),
    )
    def test(self, building_block, fg_ids, points):
        a, b, c, d = building_block.get_bonder_plane(
            fg_ids=fg_ids,
        )
        for point in points:
            product = a*point.x + b*point.y + c*point.z
            assert abs(point.distance - abs(product-d)) < 1e-16


class TestGetBonderPlaneNormal:
    def case1():
        building_block = stk.BuildingBlock(
            smiles='BrCC(CBr)(CBr)CBr',
            functional_groups=['bromine'],
        )
        bonder_ids = list(building_block.get_bonder_ids())
        coords = building_block.get_position_matrix()
        coords[bonder_ids[0]] = [1, 1, 0]
        coords[bonder_ids[1]] = [0, 0, 0.5]
        coords[bonder_ids[2]] = [0, 0, -0.5]
        coords[bonder_ids[3]] = [1, -1, 0]
        building_block.set_position_matrix(coords)
        return building_block, None, [0, 0, -1]

    @pytest.mark.parametrize(
        argnames=(
            'building_block',
            'fg_ids',
            'expected_bonder_plane_normal',
        ),
        argvalues=(
            case1(),
        ),
    )
    def test(
        self,
        building_block,
        fg_ids,
        expected_bonder_plane_normal,
    ):
        bonder_plane_normal = building_block.get_bonder_plane_normal(
            fg_ids=fg_ids,
        )
        assert np.allclose(
            a=bonder_plane_normal,
            b=expected_bonder_plane_normal,
            atol=1e-32,
        )


class TestBonderDistances:
    def case1():
        building_block = stk.BuildingBlock(
            smiles='BrCC(CBr)(CBr)CBr',
            functional_groups=['bromine'],
        )
        coords = building_block.get_position_matrix()
        bonders = enumerate(building_block.get_bonder_ids())
        for fg, bonder_id in bonders:
            coords[bonder_id] = [fg, 0, 0]
        building_block.set_position_matrix(coords)
        distances = {
            (0, 1): 1,
            (0, 2): 2,
            (0, 3): 3,
            (1, 2): 1,
            (1, 3): 2,
            (2, 3): 1,
        }
        return building_block, None, distances

    @pytest.mark.parametrize(
        argnames=(
            'building_block',
            'fg_ids',
            'expected_distances',
        ),
        argvalues=(
            case1(),
        ),
    )
    def test(
        self,
        building_block,
        fg_ids,
        expected_distances,
    ):
        distances = building_block.get_bonder_distances(
            fg_ids=fg_ids,
        )
        for fg1, fg2, distance in distances:
            assert (
                abs(expected_distances[(fg1, fg2)] - distance) < 1e-32
            )


class TestGetBonderDirectionVectors:
    def case1():
        building_block = stk.BuildingBlock(
            smiles='BrCC(CBr)(CBr)CBr',
            functional_groups=['bromine'],
        )
        coords = building_block.get_position_matrix()
        bonders = enumerate(building_block.get_bonder_ids())
        for fg_id, bonder_id in bonders:
            coords[bonder_id] = [fg_id, fg_id, fg_id]
        building_block.set_position_matrix(coords)
        expected_vectors = {
            (1, 0): [-1, -1, -1],
            (2, 0): [-2, -2, -2],
            (3, 0): [-3, -3, -3],
            (2, 1): [-1, -1, -1],
            (3, 1): [-2, -2, -2],
            (3, 2): [-1, -1, -1],
        }
        return building_block, None, expected_vectors

    @pytest.mark.parametrize(
        argnames=(
            'building_block',
            'fg_ids',
            'expected_vectors',
        ),
        argvalues=(
            case1(),
        ),
    )
    def test(
        self,
        building_block,
        fg_ids,
        expected_vectors,
    ):
        direction_vectors = (
            building_block.get_bonder_direction_vectors(fg_ids=fg_ids)
        )
        for fg1, fg2, direction in direction_vectors:
            assert np.allclose(
                a=direction,
                b=expected_vectors[(fg1, fg2)],
                atol=1e-32,
            )


class TestGetCentroidCentroidDirectionVector:
    def case1():
        building_block = stk.BuildingBlock(
            smiles='BrCC(CBr)(CBr)CBr',
            functional_groups=['bromine'],
        )
        coords = np.zeros((len(building_block.atoms), 3))
        for bonder_id in building_block.get_bonder_ids():
            coords[bonder_id] = [1, 0, 0]
        building_block.set_position_matrix(coords)
        return building_block, None, [-1, 0, 0]

    @pytest.mark.parametrize(
        argnames=(
            'building_block',
            'fg_ids',
            'expected_vector',
        ),
        argvalues=(
            case1(),
        )
    )
    def test(
        self,
        building_block,
        fg_ids,
        expected_vector,
    ):
        direction = (
            building_block.get_centroid_centroid_direction_vector(
                fg_ids=fg_ids,
            )
        )
        assert np.allclose(
            a=stk.normalize_vector(direction),
            b=expected_vector,
            atol=1e-32,
        )


def test_dump_and_load(tmpdir, building_block):
    path = str(tmpdir / 'building_block.json')
    building_block.dump(path)
    new_building_block = stk.BuildingBlock.load(path)
    assert is_equivalent_building_block(
        building_block1=building_block,
        building_block2=new_building_block,
    )
