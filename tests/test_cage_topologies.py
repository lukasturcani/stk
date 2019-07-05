import stk
import os
from os.path import join
import numpy as np


test_dir = 'cage_topology_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


def test_alignments(amine2, amine2_alt3, aldehyde3, aldehyde3_alt3):
    # Regular cage topology.
    bb_positions = {
        amine2: [0, 1, 2, 3, 4],
        amine2_alt3: [5],
        aldehyde3: [0, 1, 2],
        aldehyde3_alt3: [3]
    }

    bbs = [amine2, amine2_alt3, aldehyde3, aldehyde3_alt3]
    for fg in range(3):
        top = stk.FourPlusSix(bb_positions=bb_positions,
                              A_alignments=[0, 0, 0, fg],
                              B_alignments=[1, 1, 1, 1, 1, 1])
        c = stk.Cage(bbs, top)
        c.write(join(test_dir, f'4p6_valignment_{fg}.mol'))

    top = stk.FourPlusSix(bb_positions=bb_positions,
                          A_alignments=[0, 0, 0, 0],
                          B_alignments=[1, 1, 1, 1, 1, -1])
    c = stk.Cage(bbs, top)
    c.write(join(test_dir, f'4p6_edge_alignment.mol'))

    # No linker topology.
    bbs = [aldehyde3, aldehyde3_alt3]
    bb_positions = {
        aldehyde3: [1, 2, 3],
        aldehyde3_alt3: [0]
    }
    for fg in range(3):
        top = stk.TwoPlusTwo(bb_positions=bb_positions,
                             alignments=[fg, 0, 0, 0])
        c = stk.Cage(bbs, top)
        c.write(join(test_dir, f'2p2_valignment_{fg}.mol'))


def test_SixPlusEight(aldehyde3, amine4):
    top = stk.SixPlusEight()
    amine_fg_count = 4
    amine_count = 6
    aldehyde_count = 8

    c = stk.Cage([aldehyde3, amine4], top)
    c.write(join(test_dir, 'SixPlusEight.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    assert (c.mol.GetNumAtoms() ==
            amine4.mol.GetNumAtoms()*amine_count +
            aldehyde3.mol.GetNumAtoms()*aldehyde_count -
            c.bonds_made*3)
    assert (c.mol.GetNumBonds() ==
            amine4.mol.GetNumBonds()*amine_count +
            aldehyde3.mol.GetNumBonds()*aldehyde_count -
            c.bonds_made*2)
    assert c.topology == top
    assert c.bb_counter[amine4] == amine_count
    assert c.bb_counter[aldehyde3] == aldehyde_count


def test_TwoPlusFour(aldehyde2, amine4):
    top = stk.TwoPlusFour()
    amine_fg_count = 4
    amine_count = 2
    aldehyde_count = 4

    c = stk.Cage([aldehyde2, amine4], top)
    c.write(join(test_dir, 'TwoPlusFour.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    assert (c.mol.GetNumAtoms() ==
            amine4.mol.GetNumAtoms()*amine_count +
            aldehyde2.mol.GetNumAtoms()*aldehyde_count -
            c.bonds_made*3)
    assert (c.mol.GetNumBonds() ==
            amine4.mol.GetNumBonds()*amine_count +
            aldehyde2.mol.GetNumBonds()*aldehyde_count -
            c.bonds_made*2)
    assert c.topology == top
    assert c.bb_counter[amine4] == amine_count
    assert c.bb_counter[aldehyde2] == aldehyde_count


def test_ThreePlusSix(aldehyde2, amine4):
    top = stk.ThreePlusSix()
    amine_fg_count = 4
    amine_count = 3
    aldehyde_count = 6

    c = stk.Cage([aldehyde2, amine4], top)
    c.write(join(test_dir, 'ThreePlusSix.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    assert (c.mol.GetNumAtoms() ==
            amine4.mol.GetNumAtoms()*amine_count +
            aldehyde2.mol.GetNumAtoms()*aldehyde_count -
            c.bonds_made*3)
    assert (c.mol.GetNumBonds() ==
            amine4.mol.GetNumBonds()*amine_count +
            aldehyde2.mol.GetNumBonds()*aldehyde_count -
            c.bonds_made*2)
    assert c.topology == top
    assert c.bb_counter[amine4] == amine_count
    assert c.bb_counter[aldehyde2] == aldehyde_count


def test_FourPlusEight(aldehyde2, amine4):
    top = stk.FourPlusEight()
    amine_fg_count = 4
    amine_count = 4
    aldehyde_count = 8

    c = stk.Cage([aldehyde2, amine4], top)
    c.write(join(test_dir, 'FourPlusEight.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    assert (c.mol.GetNumAtoms() ==
            amine4.mol.GetNumAtoms()*amine_count +
            aldehyde2.mol.GetNumAtoms()*aldehyde_count -
            c.bonds_made*3)
    assert (c.mol.GetNumBonds() ==
            amine4.mol.GetNumBonds()*amine_count +
            aldehyde2.mol.GetNumBonds()*aldehyde_count -
            c.bonds_made*2)
    assert c.topology == top
    assert c.bb_counter[amine4] == amine_count
    assert c.bb_counter[aldehyde2] == aldehyde_count


def test_FivePlusTen(aldehyde2, amine4):
    top = stk.FivePlusTen()
    amine_fg_count = 4
    amine_count = 5
    aldehyde_count = 10

    c = stk.Cage([aldehyde2, amine4], top)
    c.write(join(test_dir, 'FivePlusTen.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    assert (c.mol.GetNumAtoms() ==
            amine4.mol.GetNumAtoms()*amine_count +
            aldehyde2.mol.GetNumAtoms()*aldehyde_count -
            c.bonds_made*3)
    assert (c.mol.GetNumBonds() ==
            amine4.mol.GetNumBonds()*amine_count +
            aldehyde2.mol.GetNumBonds()*aldehyde_count -
            c.bonds_made*2)
    assert c.topology == top
    assert c.bb_counter[amine4] == amine_count
    assert c.bb_counter[aldehyde2] == aldehyde_count


def test_SixPlusTwelve(aldehyde2, amine4):
    top = stk.SixPlusTwelve()
    amine_fg_count = 4
    amine_count = 6
    aldehyde_count = 12

    c = stk.Cage([aldehyde2, amine4], top)
    c.write(join(test_dir, 'SixPlusTwelve.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    assert (c.mol.GetNumAtoms() ==
            amine4.mol.GetNumAtoms()*amine_count +
            aldehyde2.mol.GetNumAtoms()*aldehyde_count -
            c.bonds_made*3)
    assert (c.mol.GetNumBonds() ==
            amine4.mol.GetNumBonds()*amine_count +
            aldehyde2.mol.GetNumBonds()*aldehyde_count -
            c.bonds_made*2)
    assert c.topology == top
    assert c.bb_counter[amine4] == amine_count
    assert c.bb_counter[aldehyde2] == aldehyde_count


def test_EightPlusSixteen(aldehyde2, amine4):
    top = stk.EightPlusSixteen()
    amine_fg_count = 4
    amine_count = 8
    aldehyde_count = 16

    c = stk.Cage([aldehyde2, amine4], top)
    c.write(join(test_dir, 'EightPlusSixteen.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    assert (c.mol.GetNumAtoms() ==
            amine4.mol.GetNumAtoms()*amine_count +
            aldehyde2.mol.GetNumAtoms()*aldehyde_count -
            c.bonds_made*3)
    assert (c.mol.GetNumBonds() ==
            amine4.mol.GetNumBonds()*amine_count +
            aldehyde2.mol.GetNumBonds()*aldehyde_count -
            c.bonds_made*2)
    assert c.topology == top
    assert c.bb_counter[amine4] == amine_count
    assert c.bb_counter[aldehyde2] == aldehyde_count


def test_TenPlusTwenty(aldehyde2, amine4):
    top = stk.TenPlusTwenty()
    amine_fg_count = 4
    amine_count = 10
    aldehyde_count = 20

    c = stk.Cage([amine4, aldehyde2], top)
    c.write(join(test_dir, 'TenPlusTwenty.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    assert (c.mol.GetNumAtoms() ==
            amine4.mol.GetNumAtoms()*amine_count +
            aldehyde2.mol.GetNumAtoms()*aldehyde_count -
            c.bonds_made*3)
    assert (c.mol.GetNumBonds() ==
            amine4.mol.GetNumBonds()*amine_count +
            aldehyde2.mol.GetNumBonds()*aldehyde_count -
            c.bonds_made*2)
    assert c.topology == top
    assert c.bb_counter[amine4] == amine_count
    assert c.bb_counter[aldehyde2] == aldehyde_count


def test_OnePlusOne(amine3, aldehyde3):
    top = stk.OnePlusOne(bb_positions={
                            0: [0],
                            1: [1]
    })
    amine_fg_count = 3
    amine_count = 1
    aldehyde_count = 1

    c = stk.Cage([aldehyde3, amine3], top)
    c.write(join(test_dir, 'OnePlusOne.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    assert (c.mol.GetNumAtoms() ==
            amine3.mol.GetNumAtoms()*amine_count +
            aldehyde3.mol.GetNumAtoms()*aldehyde_count -
            c.bonds_made*3)
    assert (c.mol.GetNumBonds() ==
            amine3.mol.GetNumBonds()*amine_count +
            aldehyde3.mol.GetNumBonds()*aldehyde_count -
            c.bonds_made*2)
    assert c.topology == top
    assert c.bb_counter[amine3] == amine_count
    assert c.bb_counter[aldehyde3] == aldehyde_count


def test_TwoPlusTwo(amine3, aldehyde3):
    top = stk.TwoPlusTwo(bb_positions={
                            0: [0, 1],
                            1: [2, 3]
    })
    amine_fg_count = 3
    amine_count = 2
    aldehyde_count = 2

    c = stk.Cage([aldehyde3, amine3], top)
    c.write(join(test_dir, 'TwoPlusTwo.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    assert (c.mol.GetNumAtoms() ==
            amine3.mol.GetNumAtoms()*amine_count +
            aldehyde3.mol.GetNumAtoms()*aldehyde_count -
            c.bonds_made*3)
    assert (c.mol.GetNumBonds() ==
            amine3.mol.GetNumBonds()*amine_count +
            aldehyde3.mol.GetNumBonds()*aldehyde_count -
            c.bonds_made*2)
    assert c.topology == top
    assert c.bb_counter[amine3] == amine_count
    assert c.bb_counter[aldehyde3] == aldehyde_count


def test_FourPlusFour(amine3, aldehyde3):
    top = stk.FourPlusFour(bb_positions={
                                0: [4, 1, 6, 3],
                                1: [0, 5, 2, 7]
    })
    amine_fg_count = 3
    amine_count = 4
    aldehyde_count = 4

    c = stk.Cage([aldehyde3, amine3], top)
    c.write(join(test_dir, 'FourPlusFour.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    assert (c.mol.GetNumAtoms() ==
            amine3.mol.GetNumAtoms()*amine_count +
            aldehyde3.mol.GetNumAtoms()*aldehyde_count -
            c.bonds_made*3)
    assert (c.mol.GetNumBonds() ==
            amine3.mol.GetNumBonds()*amine_count +
            aldehyde3.mol.GetNumBonds()*aldehyde_count -
            c.bonds_made*2)
    assert c.topology == top
    assert c.bb_counter[amine3] == amine_count
    assert c.bb_counter[aldehyde3] == aldehyde_count


def test_TwoPlusThree(amine2, aldehyde3):
    top = stk.TwoPlusThree()
    amine_fg_count = 2
    amine_count = 3
    aldehyde_count = 2

    c = stk.Cage([amine2, aldehyde3], top)
    c.write(join(test_dir, 'TwoPlusThree.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    assert (c.mol.GetNumAtoms() ==
            amine2.mol.GetNumAtoms()*amine_count +
            aldehyde3.mol.GetNumAtoms()*aldehyde_count -
            c.bonds_made*3)
    assert (c.mol.GetNumBonds() ==
            amine2.mol.GetNumBonds()*amine_count +
            aldehyde3.mol.GetNumBonds()*aldehyde_count -
            c.bonds_made*2)
    assert c.topology == top
    assert c.bb_counter[amine2] == amine_count
    assert c.bb_counter[aldehyde3] == aldehyde_count


def test_FourPlusSix(amine2, aldehyde3):
    top = stk.FourPlusSix()
    amine_fg_count = 2
    amine_count = 6
    aldehyde_count = 4

    c = stk.Cage([amine2, aldehyde3], top)
    c.write(join(test_dir, 'FourPlusSix.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    assert (c.mol.GetNumAtoms() ==
            amine2.mol.GetNumAtoms()*amine_count +
            aldehyde3.mol.GetNumAtoms()*aldehyde_count -
            c.bonds_made*3)
    assert (c.mol.GetNumBonds() ==
            amine2.mol.GetNumBonds()*amine_count +
            aldehyde3.mol.GetNumBonds()*aldehyde_count -
            c.bonds_made*2)
    assert c.topology == top
    assert c.bb_counter[amine2] == amine_count
    assert c.bb_counter[aldehyde3] == aldehyde_count


def test_multiFourPlusSix(amine2, amine2_alt1, amine2_alt2,
                          aldehyde3, aldehyde3_alt1, aldehyde3_alt2):
    top = stk.FourPlusSix(
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

    bbs = [amine2, amine2_alt1, amine2_alt2,
           aldehyde3, aldehyde3_alt1, aldehyde3_alt2]

    c = stk.Cage(bbs, top)
    c.write(join(test_dir, 'multi_FourPlusSix.mol'))

    assert c.bonds_made == amine_fg_count*6
    assert c.topology == top
    assert c.bb_counter[amine2] == 2
    assert c.bb_counter[amine2_alt1] == 3
    assert c.bb_counter[amine2_alt2] == 1
    assert c.bb_counter[aldehyde3] == 1
    assert c.bb_counter[aldehyde3_alt1] == 2
    assert c.bb_counter[aldehyde3_alt2] == 1

    assert (c.mol.GetNumAtoms() ==
            amine2.mol.GetNumAtoms()*2 +
            amine2_alt1.mol.GetNumAtoms()*3 +
            amine2_alt2.mol.GetNumAtoms()*1 +
            aldehyde3.mol.GetNumAtoms()*1 +
            aldehyde3_alt1.mol.GetNumAtoms()*2 +
            aldehyde3_alt2.mol.GetNumAtoms()*1 -
            c.bonds_made*3)
    assert (c.mol.GetNumBonds() ==
            amine2.mol.GetNumBonds()*2 +
            amine2_alt1.mol.GetNumBonds()*3 +
            amine2_alt2.mol.GetNumBonds()*1 +
            aldehyde3.mol.GetNumBonds()*1 +
            aldehyde3_alt1.mol.GetNumBonds()*2 +
            aldehyde3_alt2.mol.GetNumBonds()*1 -
            c.bonds_made*2)


def test_multiFourPlusFour(aldehyde3, aldehyde3_alt1, aldehyde3_alt2):
    top = stk.FourPlusFour(
                      bb_positions={
                          aldehyde3: [0, 1],
                          aldehyde3_alt1: [2, 3, 4, 6, 7],
                          aldehyde3_alt2: [5]
                       }
    )
    bbs = [aldehyde3, aldehyde3_alt1, aldehyde3_alt2]
    c = stk.Cage(bbs, top)
    c.write(join(test_dir, 'multi_FourPlusFour.mol'))

    assert c.bb_counter[aldehyde3] == 2
    assert c.bb_counter[aldehyde3_alt1] == 5
    assert c.bb_counter[aldehyde3_alt2] == 1


def test_FourPlusSix2(amine2, aldehyde3):
    top = stk.FourPlusSix2()
    amine_fg_count = 2
    amine_count = 6
    aldehyde_count = 4

    c = stk.Cage([amine2, aldehyde3], top)
    c.write(join(test_dir, 'FourPlusSix2.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    assert (c.mol.GetNumAtoms() ==
            amine2.mol.GetNumAtoms()*amine_count +
            aldehyde3.mol.GetNumAtoms()*aldehyde_count -
            c.bonds_made*3)
    assert (c.mol.GetNumBonds() ==
            amine2.mol.GetNumBonds()*amine_count +
            aldehyde3.mol.GetNumBonds()*aldehyde_count -
            c.bonds_made*2)
    assert c.topology == top
    assert c.bb_counter[amine2] == amine_count
    assert c.bb_counter[aldehyde3] == aldehyde_count


def test_SixPlusNine(amine2, aldehyde3):
    top = stk.SixPlusNine()
    amine_fg_count = 2
    amine_count = 9
    aldehyde_count = 6

    c = stk.Cage([amine2, aldehyde3], top)
    c.write(join(test_dir, 'SixPlusNine.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    assert (c.mol.GetNumAtoms() ==
            amine2.mol.GetNumAtoms()*amine_count +
            aldehyde3.mol.GetNumAtoms()*aldehyde_count -
            c.bonds_made*3)
    assert (c.mol.GetNumBonds() ==
            amine2.mol.GetNumBonds()*amine_count +
            aldehyde3.mol.GetNumBonds()*aldehyde_count -
            c.bonds_made*2)
    assert c.topology == top
    assert c.bb_counter[amine2] == amine_count
    assert c.bb_counter[aldehyde3] == aldehyde_count


def test_EightPlusTwelve(amine2, aldehyde3):
    top = stk.EightPlusTwelve()
    amine_fg_count = 2
    amine_count = 12
    aldehyde_count = 8

    c = stk.Cage([amine2, aldehyde3], top)
    c.write(join(test_dir, 'EightPlusTwelve.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    assert (c.mol.GetNumAtoms() ==
            amine2.mol.GetNumAtoms()*amine_count +
            aldehyde3.mol.GetNumAtoms()*aldehyde_count -
            c.bonds_made*3)
    assert (c.mol.GetNumBonds() ==
            amine2.mol.GetNumBonds()*amine_count +
            aldehyde3.mol.GetNumBonds()*aldehyde_count -
            c.bonds_made*2)
    assert c.topology == top
    assert c.bb_counter[amine2] == amine_count
    assert c.bb_counter[aldehyde3] == aldehyde_count


def test_Dodecahedron(amine2, aldehyde3):
    top = stk.Dodecahedron()
    amine_fg_count = 2
    amine_count = 30
    aldehyde_count = 20

    c = stk.Cage([amine2, aldehyde3], top)
    c.write(join(test_dir, 'Dodecahedron.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    assert (c.mol.GetNumAtoms() ==
            amine2.mol.GetNumAtoms()*amine_count +
            aldehyde3.mol.GetNumAtoms()*aldehyde_count -
            c.bonds_made*3)
    assert (c.mol.GetNumBonds() ==
            amine2.mol.GetNumBonds()*amine_count +
            aldehyde3.mol.GetNumBonds()*aldehyde_count -
            c.bonds_made*2)
    assert c.topology == top
    assert c.bb_counter[amine2] == amine_count
    assert c.bb_counter[aldehyde3] == aldehyde_count


def test_Icosahedron(amine2, aldehyde5):
    top = stk.Icosahedron()
    amine_fg_count = 2
    amine_count = 30
    aldehyde_count = 12

    c = stk.Cage([amine2, aldehyde5], top)
    c.write(join(test_dir, 'Icosahedron.mol'))

    assert c.bonds_made == amine_fg_count*amine_count
    assert (c.mol.GetNumAtoms() ==
            amine2.mol.GetNumAtoms()*amine_count +
            aldehyde5.mol.GetNumAtoms()*aldehyde_count -
            c.bonds_made*3)
    assert (c.mol.GetNumBonds() ==
            amine2.mol.GetNumBonds()*amine_count +
            aldehyde5.mol.GetNumBonds()*aldehyde_count -
            c.bonds_made*2)
    assert c.topology == top
    assert c.bb_counter[amine2] == amine_count
    assert c.bb_counter[aldehyde5] == aldehyde_count


def test_multiconformer(tmp_amine2, tmp_aldehyde3):
    top = stk.FourPlusSix()
    amine_fg_count = 2
    amine_count = 6
    aldehyde_count = 4

    c = stk.Cage([tmp_amine2, tmp_aldehyde3],
                 stk.FourPlusSix(),
                 bb_conformers=[0, 0])
    c.add_conformer([1, 0])
    c.add_conformer([0, 1])
    c.add_conformer([1, 1])

    c.write(join(test_dir, 'FourPlusSix_conf1.mol'), conformer=0)
    c.write(join(test_dir, 'FourPlusSix_conf2.mol'), conformer=1)
    c.write(join(test_dir, 'FourPlusSix_conf3.mol'), conformer=2)
    c.write(join(test_dir, 'FourPlusSix_conf4.mol'), conformer=3)

    assert c.bonds_made == amine_fg_count*amine_count
    assert (c.mol.GetNumAtoms() ==
            tmp_amine2.mol.GetNumAtoms()*amine_count +
            tmp_aldehyde3.mol.GetNumAtoms()*aldehyde_count -
            c.bonds_made*3)
    assert (c.mol.GetNumBonds() ==
            tmp_amine2.mol.GetNumBonds()*amine_count +
            tmp_aldehyde3.mol.GetNumBonds()*aldehyde_count -
            c.bonds_made*2)
    assert c.topology == top
    assert c.bb_counter[tmp_amine2] == amine_count
    assert c.bb_counter[tmp_aldehyde3] == aldehyde_count


def test_cage_complex(amine2, amine2_alt1, aldehyde3, chained_c60):
    c = stk.Cage([amine2, amine2_alt1, aldehyde3],
                 stk.FourPlusSix(bb_positions={
                     amine2: [5],
                     amine2_alt1: [0, 1, 2, 3, 4],
                     aldehyde3: [0, 1, 2, 3]
                 }))

    n = 3
    for i in range(n):
        cage_complex = stk.CageComplex(
            [c, chained_c60],
            stk.CageWithGuest(axis=[1, 0, 0],
                              angle=2*np.pi*i/n,
                              displacement=[2*i, 0, 0])
        )
        cage_complex.write(join(test_dir, f'cage_with_guest_{i}.mol'))
