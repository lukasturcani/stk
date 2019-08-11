import stk
import os
from os.path import join
import itertools as it


test_dir = 'cage_topology_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


def _add_test_attrs(cage):
    cage.test_attr1 = 'something'
    cage.test_attr2 = 12
    cage.test_attr3 = ['12', 'something', 21]
    cage.test_attr4 = 'skip'

    bb1, bb2 = cage.building_block_vertices
    bb1.test_attr1 = 1232
    bb2.test_attr5 = 'alpha'

    # Add some custom atom properties.
    cage.atoms[0].some_prop = 'custom atom prop'
    # Add some custom bond properties
    cage.bonds[2].other_prop = 1999


def _test_func_groups(cage, loaded):
    fgs = it.zip_longest(cage.func_groups, loaded.func_groups)
    for fg1, fg2 in fgs:
        atoms = it.zip_longest(fg1.atoms, fg2.atoms)
        bonders = it.zip_longest(fg1.bonders, fg2.bonders)
        deleters = it.zip_longest(fg1.deleters, fg2.deleters)
        for a1, a2 in it.chain(atoms, bonders, deleters):
            assert a1.__class__ is a2.__class__
            assert a1.id is a1.id


def _test_atoms(cage, loaded):
    for a1, a2 in zip(cage.atoms, loaded.atoms):
        assert a1.__class__ is a2.__class__
        d1, d2 = vars(a1), vars(a2)
        bb1, bb2 = d1.pop('building_block'), d2.pop('building_block')
        assert d1 == d2
        assert bb1.is_identical(bb2)


def _test_bonds(cage, loaded):
    for b1, b2 in zip(cage.bonds, loaded.bonds):
        assert b1.__class__ is b2.__class__
        d1, d2 = vars(b1), vars(b2)
        assert repr(d1.pop('atom1')) == repr(d2.pop('atom1'))
        assert repr(d1.pop('atom2')) == repr(d2.pop('atom2'))
        assert d1 == d2


def _test_attrs(cage, loaded):
    assert cage.test_attr1 == loaded.test_attr1
    assert cage.test_attr2 == loaded.test_attr2
    assert cage.test_attr3 == loaded.test_attr3
    assert not hasattr(loaded, 'test_attr4')


def _test_bbs(cage, loaded):
    bbs1 = list(cage.building_block_vertices.keys())
    bbs2 = list(loaded.building_block_vertices.keys())
    for bb1, bb2 in it.zip_longest(bbs1, bbs2):
        assert bb1.is_identical(bb2)
        assert bb1 is not bb2
        bb1_count = cage.building_block_counter[bb1]
        bb2_count = loaded.building_block_counter[bb2]
        assert bb1_count == bb2_count

    assert bbs2[0].test_attr1 == 1232
    assert bbs2[1].test_attr5 == 'alpha'


def _test_dump_and_load(cage):
    path = join(
        test_dir, f'{cage.topology_graph.__class__.__name__}.dump'
    )
    _add_test_attrs(cage)
    cage.dump(
        path=path,
        include_attrs=[
            'test_attr1',
            'test_attr2',
            'test_attr3',
            'test_attr5'
        ],
        ignore_missing_attrs=True
    )
    loaded = stk.Molecule.load(path)

    assert cage.__class__ is loaded.__class__
    assert loaded is not cage
    _test_func_groups(cage, loaded)
    _test_atoms(cage, loaded)
    _test_bonds(cage, loaded)

    assert repr(loaded.topology_graph) == repr(cage.topology_graph)
    assert (
        len(loaded.construction_bonds) == len(cage.construction_bonds)
    )
    _test_attrs(cage, loaded)
    _test_bbs(cage, loaded)
    mol3 = stk.Molecule.load(path, use_cache=True)
    mol4 = stk.Molecule.load(path, use_cache=True)
    assert mol3 is mol4


def test_topologies(
    tmp_six_plus_eight,
    tmp_one_plus_one,
    tmp_two_plus_two,
    tmp_four_plus_four,
    tmp_twelve_plus_thirty,
    tmp_two_plus_four,
    tmp_three_plus_six,
    tmp_four_plus_eight,
    tmp_five_plus_ten,
    tmp_six_plus_twelve,
    tmp_eight_plus_sixteen,
    tmp_ten_plus_twenty,
    tmp_two_plus_three,
    tmp_four_plus_six,
    tmp_four_plus_six2,
    tmp_six_plus_nine,
    tmp_eight_plus_twelve,
    tmp_twenty_plus_thirty
):
    cages = (
        (tmp_six_plus_eight, 6, 8),
        (tmp_one_plus_one, 1, 1),
        (tmp_two_plus_two, 2, 2),
        (tmp_four_plus_four, 4, 4),
        (tmp_twelve_plus_thirty, 12, 30),
        (tmp_two_plus_four, 2, 4),
        (tmp_three_plus_six, 3, 6),
        (tmp_four_plus_eight, 4, 8),
        (tmp_five_plus_ten, 5, 10),
        (tmp_six_plus_twelve, 6, 12),
        (tmp_eight_plus_sixteen, 8, 16),
        (tmp_ten_plus_twenty, 10, 20),
        (tmp_two_plus_three, 2, 3),
        (tmp_four_plus_six, 4, 6),
        (tmp_four_plus_six2, 4, 6),
        (tmp_six_plus_nine, 6, 9),
        (tmp_eight_plus_twelve, 8, 12),
        (tmp_twenty_plus_thirty, 20, 30),
    )
    for cage, num_expected_bb1s, num_expected_bb2s in cages:
        _test_construction(cage, num_expected_bb1s, num_expected_bb2s)
        _test_dump_and_load(cage)


def test_alignments(amine2, amine2_alt3, aldehyde3, aldehyde3_alt3):

    building_blocks = [amine2, amine2_alt3, aldehyde3, aldehyde3_alt3]
    for fg in range(3):
        v4 = stk.cage.FourPlusSix.vertices[3]
        tetrahedron = stk.cage.FourPlusSix(
            vertex_alignments={
                v4: v4.edges[fg]
            },
        )

        building_block_vertices = {
            amine2: tetrahedron.vertices[4:9],
            amine2_alt3: tetrahedron.vertices[9:],
            aldehyde3: tetrahedron.vertices[:3],
            aldehyde3_alt3: tetrahedron.vertices[3:4]
        }
        c = stk.ConstructedMolecule(
            buliding_blocks=building_blocks,
            topology_graph=tetrahedron,
            building_block_vertices=building_block_vertices
        )
        c.write(join(test_dir, f'4p6_valignment_{fg}.mol'))

    v10 = stk.cage.FourPlusSix.vertices[9]
    tetrahedron = stk.FourPlusSix(
        vertex_alignments={
            v10: v10.edges[1]
        }
    )
    building_block_vertices = {
        amine2: tetrahedron.vertices[4:9],
        amine2_alt3: tetrahedron.vertices[9:],
        aldehyde3: tetrahedron.vertices[:3],
        aldehyde3_alt3: tetrahedron.vertices[3:4]
    }
    c = stk.ConstructedMolecule(
        buliding_blocks=building_blocks,
        topology_graph=tetrahedron,
        building_block_vertices=building_block_vertices
    )
    c.write(join(test_dir, f'4p6_edge_alignment.mol'))


def _test_construction(cage, num_expected_bb1s, num_expected_bb2s):
    cage.write(join(test_dir, f'{cage.__class__.__name__}.mol'))
    bb1, bb2 = sorted(
        cage.building_block_vertices,
        key=lambda bb: len(bb.func_groups),
        reverse=True
    )
    num_bb1s = cage.building_block_counter[bb1]
    num_bb2s = cage.building_block_counter[bb2]
    num_bb1_deleters = sum(len(fg.deleters) for fg in bb1.func_groups)
    num_bb2_deleters = sum(len(fg.deleters) for fg in bb2.func_groups)

    # Check that the correct number of bonds got made.
    assert (
        len(cage.construction_bonds) == len(cage.topology_graph.edges)
    )
    assert (
        len(cage.construction_bonds) == num_bb1s * len(bb1.func_groups)
    )
    assert (
        len(cage.construction_bonds) == num_bb2s * len(bb2.func_groups)
    )
    # Check correct total number of atoms.
    assert (
        len(cage.atoms) ==
        len(bb1.atoms)*num_bb1s
        + len(bb2.atoms)*num_bb2s
        - num_bb1_deleters*num_bb1s
        - num_bb2_deleters*num_bb2s
    )
    # Check correct total number of bonds.
    assert (
        len(cage.bonds) ==
        len(bb1.bonds)*num_bb1s
        + len(bb2.bonds)*num_bb2s
        + len(cage.construction_bonds)
        - num_bb1_deleters*num_bb1s
        - num_bb2_deleters*num_bb2s
    )
    # Check building block counts.
    assert num_bb1s == num_expected_bb1s
    assert num_bb2s == num_expected_bb2s
    # Check window attributes got added
    assert cage.num_windows == cage.topology_graph.num_windows
    assert (
        cage.num_window_types == cage.topology_graph.num_window_types
    )


def test_multicage(
    amine2, amine2_alt1, amine2_alt2,
    aldehyde3, aldehyde3_alt1, aldehyde3_alt2
):
    tetrahedron = stk.FourPlusSix(
                      bb_positions={
                          0: [0, 1],
                          1: [2, 3, 4],
                          2: [5],
                          3: [0],
                          4: [1, 2],
                          5: [3]
                       }
    )

    amine_fg_count = 2

    building_blocks = [
        amine2, amine2_alt1, amine2_alt2,
        aldehyde3, aldehyde3_alt1, aldehyde3_alt2
    ]

    c = stk.ConstructedMolecule(
        buliding_blocks=building_blocks,
        topology_graph=tetrahedron
    )
    c.write(join(test_dir, 'multi_FourPlusSix.mol'))

    assert c.bonds_made == amine_fg_count*6
    assert repr(c.topology) == repr(tetrahedron)
    assert c.building_block_counter[amine2] == 2
    assert c.building_block_counter[amine2_alt1] == 3
    assert c.building_block_counter[amine2_alt2] == 1
    assert c.building_block_counter[aldehyde3] == 1
    assert c.building_block_counter[aldehyde3_alt1] == 2
    assert c.building_block_counter[aldehyde3_alt2] == 1

    num_expected_atoms = (
        len(amine2.atoms)*2 +
        len(amine2_alt1.atoms)*3 +
        len(amine2_alt2.atoms) +
        len(aldehyde3.atoms) +
        len(aldehyde3_alt1.atoms)*2 +
        len(aldehyde3_alt2.atoms) -
        c.bonds_made*3
    )
    assert len(c.atoms) == num_expected_atoms
    num_expected_bonds = (
        len(amine2.bonds)*2 +
        len(amine2_alt1.bonds)*3 +
        len(amine2_alt2.bonds)*1 +
        len(aldehyde3.bonds)*1 +
        len(aldehyde3_alt1.bonds)*2 +
        len(aldehyde3_alt2.bonds)*1 -
        c.bonds_made*2
    )
    assert len(c.bonds) == num_expected_bonds
