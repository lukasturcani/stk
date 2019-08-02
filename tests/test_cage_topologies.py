import stk
import os
from os.path import join
import numpy as np


test_dir = 'cage_topology_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


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


def test_SixPlusEight(aldehyde3, amine4):
    rhombic_dodecahedron = stk.cage.SixPlusEight()
    amine_fg_count = 4
    amine_count = 6
    aldehyde_count = 8

    c = stk.ConstructedMolecule(
        buliding_blocks=[aldehyde3, amine4],
        topology_graph=rhombic_dodecahedron
    )
    c.write(join(test_dir, 'SixPlusEight.mol'))

    assert c.bonds_made == amine_fg_count * amine_count
    num_expected_atoms = (
        len(amine4.atoms)*amine_count +
        len(aldehyde3.atoms)*aldehyde_count -
        c.bonds_made*3
    )
    assert len(c.atoms) == num_expected_atoms
    num_expected_bonds = (
        len(amine4.bonds)*amine_count +
        len(aldehyde3.bonds)*aldehyde_count -
        c.bonds_made*2
    )
    assert len(c.bonds) == num_expected_bonds
    assert c.topology_graph == rhombic_dodecahedron
    assert c.building_block_counter[amine4] == amine_count
    assert c.building_block_counter[aldehyde3] == aldehyde_count


def test_TwoPlusFour(aldehyde2, amine4):
    capsule = stk.cage.TwoPlusFour()
    amine_fg_count = 4
    amine_count = 2
    aldehyde_count = 4

    c = stk.ConstructedMolecule(
        buliding_blocks=[aldehyde2, amine4],
        topology_graph=capsule
    )
    c.write(join(test_dir, 'TwoPlusFour.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    num_expected_atoms = (
        len(amine4.atoms)*amine_count +
        len(aldehyde2.atoms)*aldehyde_count -
        c.bonds_made*3
    )
    assert len(c.atoms) == num_expected_atoms
    num_expected_bonds = (
        len(amine4.bonds)*amine_count +
        len(aldehyde2.bonds)*aldehyde_count -
        c.bonds_made*2
    )
    assert len(c.bonds) == num_expected_bonds
    assert c.topology_graph == capsule
    assert c.building_block_counter[amine4] == amine_count
    assert c.building_block_counter[aldehyde2] == aldehyde_count


def test_ThreePlusSix(aldehyde2, amine4):
    top = stk.cage.ThreePlusSix()
    amine_fg_count = 4
    amine_count = 3
    aldehyde_count = 6

    c = stk.ConstructedMolecule(
        building_blocks=[aldehyde2, amine4],
        topology_graph=topology_graph
    )
    c.write(join(test_dir, 'ThreePlusSix.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    num_expected_atoms = (
        amine4.mol.GetNumAtoms()*amine_count +
        aldehyde2.mol.GetNumAtoms()*aldehyde_count -
        c.bonds_made*3
    )
    assert len(c.atoms) == num_expected_atoms
    num_expected_bonds = (
        len(amine4.bonds)*amine_count +
        len(aldehyde2.bonds)*aldehyde_count -
        c.bonds_made*2
    )
    assert len(c.bonds) == num_expected_bonds
    assert c.topology_graph == a
    assert c.building_block_counter[amine4] == amine_count
    assert c.building_block_counter[aldehyde2] == aldehyde_count


def test_FourPlusEight(aldehyde2, amine4):
    top = stk.cage.FourPlusEight()
    amine_fg_count = 4
    amine_count = 4
    aldehyde_count = 8

    c = stk.ConstructedMolecule(
        buliding_blocks=[aldehyde2, amine4],
        topology_graph=a
    )
    c.write(join(test_dir, 'FourPlusEight.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    num_expected_atoms = (
        len(amine4.atoms)*amine_count +
        len(aldehyde2.atoms)*aldehyde_count -
        c.bonds_made*3
    )
    assert len(c.atoms) == num_expected_atoms
    num_expected_bonds = (
        len(amine4.bonds)*amine_count +
        len(aldehyde2.bonds)*aldehyde_count -
        c.bonds_made*2
    )
    assert len(c.bonds) == num_expected_bonds
    assert c.topology_graph == a
    assert c.building_block_counter[amine4] == amine_count
    assert c.building_block_counter[aldehyde2] == aldehyde_count


def test_FivePlusTen(aldehyde2, amine4):
    top = stk.cage.FivePlusTen()
    amine_fg_count = 4
    amine_count = 5
    aldehyde_count = 10

    c = stk.ConstructedMolecule(
        buliding_blocks=[aldehyde2, amine4],
        topology_graph=a
    )

    c.write(join(test_dir, 'FivePlusTen.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    num_expected_atoms = (
        len(amine4.atoms)*amine_count +
        len(aldehyde2.atoms)*aldehyde_count -
        c.bonds_made*3
    )
    assert len(c.atoms) == num_expected_atoms
    num_expected_bonds = (
        len(amine4.bonds)*amine_count +
        len(aldehyde2.bonds)*aldehyde_count -
        c.bonds_made*2
    )
    assert len(c.bonds) == num_expected_bonds
    assert c.topology_graph == a
    assert c.building_block_counter[amine4] == amine_count
    assert c.building_block_counter[aldehyde2] == aldehyde_count


def test_SixPlusTwelve(aldehyde2, amine4):
    top = stk.cage.SixPlusTwelve()
    amine_fg_count = 4
    amine_count = 6
    aldehyde_count = 12

    c = stk.ConstructedMolecule(
        buliding_blocks=[aldehyde2, amine4],
        topology_graph=a)

    c.write(join(test_dir, 'SixPlusTwelve.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    num_expected_atoms = (
        len(amine4.atoms)*amine_count +
        len(aldehyde2.atoms)*aldehyde_count -
        c.bonds_made*3
    )
    assert len(c.atoms) == num_expected_atoms
    num_expected_bonds = (
        len(amine4.bonds)*amine_count +
        len(aldehyde2.bonds)*aldehyde_count -
        c.bonds_made*2
    )
    assert len(c.bonds) == num_expected_bonds
    assert c.topology_graph == a
    assert c.building_block_counter[amine4] == amine_count
    assert c.building_block_counter[aldehyde2] == aldehyde_count


def test_EightPlusSixteen(aldehyde2, amine4):
    top = stk.cage.EightPlusSixteen()
    amine_fg_count = 4
    amine_count = 8
    aldehyde_count = 16

    c = stk.ConstructedMolecule(
        buliding_blocks=[aldehyde2, amine4],
        topology_graph=a
    )

    c.write(join(test_dir, 'EightPlusSixteen.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    num_expected_atoms = (
        len(amine4.atoms)*amine_count +
        len(aldehyde2.atoms)*aldehyde_count -
        c.bonds_made*3
    )
    assert len(c.atoms) == num_expected_atoms
    num_expected_bonds = (
        len(amine4.bonds)*amine_count +
        len(aldehyde2.bonds)*aldehyde_count -
        c.bonds_made*2
    )
    assert len(c.bonds) == num_expected_bonds
    assert c.topology == a
    assert c.building_block_counter[amine4] == amine_count
    assert c.building_block_counter[aldehyde2] == aldehyde_count


def test_TenPlusTwenty(aldehyde2, amine4):
    top = stk.cage.TenPlusTwenty()
    amine_fg_count = 4
    amine_count = 10
    aldehyde_count = 20

    c = stk.ConstructedMolecule(
        buliding_blocks=[amine4, aldehyde2],
        topology_graph=a
    )
    c.write(join(test_dir, 'TenPlusTwenty.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    assert len(c.atoms) == (
        len(amine4.atoms)*amine_count +
        len(aldehyde2.atoms)*aldehyde_count -
        c.bonds_made*3
    )
    num_expected_bonds = (
        len(amine4.bonds)*amine_count +
        len(aldehyde2.bonds)*aldehyde_count -
        c.bonds_made*2
    )
    assert len(c.bonds) == num_expected_bonds
    assert c.topology_graph == a
    assert c.building_block_counter[amine4] == amine_count
    assert c.building_block_counter[aldehyde2] == aldehyde_count


def test_OnePlusOne(amine3, aldehyde3):
    capsule = stk.cage.OnePlusOne(bb_positions={
                            0: [0],
                            1: [1]
    })
    amine_fg_count = 3
    amine_count = 1
    aldehyde_count = 1

    c = stk.Cage([aldehyde3, amine3], capsule)
    c.write(join(test_dir, 'OnePlusOne.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    num_expected_atoms = (
        len(amine3.atoms)*amine_count +
        len(aldehyde3.atoms)*aldehyde_count -
        c.bonds_made*3
    )
    assert len(c.atoms) == num_expected_atoms
    num_expected_bonds = (
        len(amine3.bonds)*amine_count +
        len(aldehyde3.bonds)*aldehyde_count -
        c.bonds_made*2
    )
    assert len(c.bonds) == num_expected_bonds
    assert c.topology_graph == capsule
    assert c.building_block_counter[amine3] == amine_count
    assert c.building_block_counter[aldehyde3] == aldehyde_count


def test_TwoPlusTwo(amine3, aldehyde3):
    tetrahedron = stk.cage.TwoPlusTwo(bb_positions={
                            0: [0, 1],
                            1: [2, 3]
    })
    amine_fg_count = 3
    amine_count = 2
    aldehyde_count = 2

    c = stk.ConstructedMolecule(
        building_blocks=[aldehyde3, amine3],
        topology_graph=tetrahedron
    )
    c.write(join(test_dir, 'TwoPlusTwo.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    num_expected_atoms = (
        len(amine3.atoms)*amine_count +
        len(aldehyde3.atoms)*aldehyde_count -
        c.bonds_made*3
    )
    assert len(c.atoms) == num_expected_atoms
    num_expected_bonds = (
        len(amine3.bonds)*amine_count +
        len(aldehyde3.bonds)*aldehyde_count -
        c.bonds_made*2
    )
    assert len(c.bonds) == num_expected_bonds
    assert repr(c.topology_graph) == repr(tetrahedron)
    assert c.building_block_counter[amine3] == amine_count
    assert c.building_block_counter[aldehyde3] == aldehyde_count


def test_FourPlusFour(amine3, aldehyde3):
    cube = stk.cage.FourPlusFour(bb_positions={
                                0: [4, 1, 6, 3],
                                1: [0, 5, 2, 7]
    })
    amine_fg_count = 3
    amine_count = 4
    aldehyde_count = 4

    c = stk.ConstructedMolecule(
        building_blocks=[aldehyde3, amine3],
        topology_graph=cube
    )
    c.write(join(test_dir, 'FourPlusFour.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    num_expected_atoms = (
        len(amine3.atoms)*amine_count +
        len(aldehyde3.atoms)*aldehyde_count -
        c.bonds_made*3
    )
    assert len(c.atoms) == num_expected_atoms
    num_expected_bonds = (
        len(amine3.bonds)*amine_count +
        len(aldehyde3.bonds)*aldehyde_count -
        c.bonds_made*2
    )
    assert len(c.bonds) == num_expected_bonds
    assert repr(c.topology_graph) == repr(cube)
    assert c.building_block_counter[amine3] == amine_count
    assert c.building_block_counter[aldehyde3] == aldehyde_count


def test_TwoPlusThree(amine2, aldehyde3):
    capsule = stk.cage.TwoPlusThree()
    amine_fg_count = 2
    amine_count = 3
    aldehyde_count = 2

    c = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde3],
        topology_graph=capsule
    )
    c.write(join(test_dir, 'TwoPlusThree.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    num_expected_atoms = (
        len(amine2.atoms)*amine_count +
        len(aldehyde3.atoms)*aldehyde_count -
        c.bonds_made*3
    )
    assert len(c.atoms) == num_expected_atoms
    num_expected_bonds = (
        len(amine2.bonds)*amine_count +
        len(aldehyde3.bonds)*aldehyde_count -
        c.bonds_made*2
    )
    assert len(c.bonds) == num_expected_bonds
    assert c.topology_graph == capsule
    assert c.building_block_counter[amine2] == amine_count
    assert c.building_block_counter[aldehyde3] == aldehyde_count


def test_FourPlusSix(amine2, aldehyde3):
    tetrahedron = stk.cage.FourPlusSix()
    amine_fg_count = 2
    amine_count = 6
    aldehyde_count = 4

    c = stk.ConstructedMolecule(
        buliding_blocks=[amine2, aldehyde3],
        topology_graph=tetrahedron
    )
    c.write(join(test_dir, 'FourPlusSix.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    num_expected_atoms = (
        len(amine2.atoms)*amine_count +
        len(aldehyde3.atoms)*aldehyde_count -
        c.bonds_made*3
    )
    assert len(c.atoms) == num_expected_atoms
    num_expected_bonds = (
        len(amine2.bonds)*amine_count +
        len(aldehyde3.bonds)*aldehyde_count -
        c.bonds_made*2
    )
    assert len(c.bonds) == num_expected_bonds
    assert repr(c.topology_graph) == repr(tetrahedron)
    assert c.building_block_counter[amine2] == amine_count
    assert c.building_block_counter[aldehyde3] == aldehyde_count


def test_multiFourPlusSix(amine2, amine2_alt1, amine2_alt2,
                          aldehyde3, aldehyde3_alt1, aldehyde3_alt2):
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


def test_multiFourPlusFour(aldehyde3, aldehyde3_alt1, aldehyde3_alt2):
    tetrahedron = stk.cage.FourPlusFour(
                      bb_positions={
                          aldehyde3: [0, 1],
                          aldehyde3_alt1: [2, 3, 4, 6, 7],
                          aldehyde3_alt2: [5]
                       }
    )
    building_blocks = [aldehyde3, aldehyde3_alt1, aldehyde3_alt2]
    c = stk.ConstructedMolecule(
        building_blocks=building_blocks,
        tetrahedron=tetrahedron
    )
    c.write(join(test_dir, 'multi_FourPlusFour.mol'))

    assert c.building_block_counter[aldehyde3] == 2
    assert c.building_block_counter[aldehyde3_alt1] == 5
    assert c.building_block_counter[aldehyde3_alt2] == 1


def test_FourPlusSix2(amine2, aldehyde3):
    top = stk.cage.FourPlusSix2()
    amine_fg_count = 2
    amine_count = 6
    aldehyde_count = 4

    c = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde3],
        topology_graph=a
    )
    c.write(join(test_dir, 'FourPlusSix2.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    num_expected_atoms = (
        len(amine2.atoms)*amine_count +
        len(aldehyde3.atoms)*aldehyde_count -
        c.bonds_made*3
    )
    assert len(c.atoms) == num_expected_atoms
    num_expected_bonds = (
        len(amine2.bonds)*amine_count +
        len(aldehyde3.bonds)*aldehyde_count -
        c.bonds_made*2
    )
    assert len(c.bonds) == num_expected_bonds
    assert repr(c.topology_graph) == repr(a)
    assert c.building_block_counter[amine2] == amine_count
    assert c.building_block_counter[aldehyde3] == aldehyde_count


def test_SixPlusNine(amine2, aldehyde3):
    triangular_prism = stk.cage.SixPlusNine()
    amine_fg_count = 2
    amine_count = 9
    aldehyde_count = 6

    c = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde3],
        topology_graph=triangular_prism
    )
    c.write(join(test_dir, 'SixPlusNine.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    num_expected_atoms = (
        len(amine2.atoms)*amine_count +
        len(aldehyde3.atoms)*aldehyde_count -
        c.bonds_made*3
    )
    assert len(c.atoms) == num_expected_atoms
    num_expected_bonds = (
        len(amine2.bonds)*amine_count +
        len(aldehyde3.bonds)*aldehyde_count -
        c.bonds_made*2
    )
    assert len(c.bonds) == num_expected_bonds
    assert c.topology == triangular_prism
    assert c.building_block_counter[amine2] == amine_count
    assert c.building_block_counter[aldehyde3] == aldehyde_count


def test_EightPlusTwelve(amine2, aldehyde3):
    cube = stk.cage.EightPlusTwelve()
    amine_fg_count = 2
    amine_count = 12
    aldehyde_count = 8

    c = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde3],
        topology_graph=cube
    )
    c.write(join(test_dir, 'EightPlusTwelve.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    num_expected_atoms = (
        len(amine2.atoms)*amine_count +
        len(aldehyde3.atoms)*aldehyde_count -
        c.bonds_made*3
    )
    assert len(c.atoms) == num_expected_atoms
    num_expected_bonds = (
        len(amine2.bonds)*amine_count +
        len(aldehyde3.bonds)*aldehyde_count -
        c.bonds_made*2
    )
    assert len(c.bonds) == num_expected_bonds
    assert repr(c.topology_graph) == repr(cube)
    assert c.building_block_counter[amine2] == amine_count
    assert c.building_block_counter[aldehyde3] == aldehyde_count


def test_TwentyPlusThirty(amine2, aldehyde3):
    dodecahedron = stk.cage.TenPlusTwenty()
    amine_fg_count = 2
    amine_count = 30
    aldehyde_count = 20

    c = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde3],
        topology_graph=dodecahedron
    )
    c.write(join(test_dir, 'Dodecahedron.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    num_expected_atoms = (
        len(amine2.atoms)*amine_count +
        len(aldehyde3.atoms)*aldehyde_count -
        c.bonds_made*3
    )
    assert len(c.atoms) == num_expected_atoms
    num_expected_bonds = (
        len(amine2.bonds)*amine_count +
        len(aldehyde3.bonds)*aldehyde_count -
        c.bonds_made*2
    )
    assert len(c.bonds) == num_expected_bonds
    assert repr(c.topology_graph) == repr(dodecahedron)
    assert c.building_block_counter[amine2] == amine_count
    assert c.building_block_counter[aldehyde3] == aldehyde_count


def test_A(amine2, aldehyde5):
    icosahedron = stk.cage.A()
    amine_fg_count = 2
    amine_count = 30
    aldehyde_count = 12

    c = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde5],
        topology_graph=icosahedron
    )
    c.write(join(test_dir, 'Icosahedron.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    num_expected_atoms = (
        len(amine2.atoms)*amine_count +
        len(aldehyde5.atoms)*aldehyde_count -
        c.bonds_made*3
    )
    assert len(c.atoms) == num_expected_atoms
    num_expected_bonds = (
        len(amine2.bonds)*amine_count +
        len(aldehyde5.bonds)*aldehyde_count -
        c.bonds_made*2
    )
    assert len(c.bonds) == num_expected_bonds
    assert repr(c.topology_graph) == repr(icosahedron)
    assert c.building_block_counter[amine2] == amine_count
    assert c.building_block_counter[aldehyde5] == aldehyde_count
