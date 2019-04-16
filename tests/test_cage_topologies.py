import stk
import os
from os.path import join


test_dir = 'cage_topology_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


def test_bb_positions():
    assert False


def test_alignements():
    assert False


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


def test_multiconformer(tmp_amine2, tmp_aldehyde3):
    top = stk.FourPlusSix()
    amine_fg_count = 2
    amine_count = 6
    aldehyde_count = 4

    # Add conformers.
    tmp_amine2.mol.AddConformer(
                        tmp_amine2.mol.GetConformer(),
                        True)
    tmp_aldehyde3.mol.AddConformer(
                        tmp_aldehyde3.mol.GetConformer(),
                        True)

    # Give conformers distinct geometries.
    tmp_amine2.set_position_from_matrix(
        pos_mat=tmp_amine2.mol.GetConformer().GetPositions().T*4,
        conformer=1)
    tmp_aldehyde3.set_position_from_matrix(
        pos_mat=tmp_aldehyde3.mol.GetConformer().GetPositions().T*4,
        conformer=1)

    c = stk.Cage([tmp_amine2, tmp_aldehyde3],
                 stk.FourPlusSix(),
                 bb_conformers=[0, 0])
    c.add_conformer([1, 0])
    c.add_conformer([0, 1])
    c.add_conformer([1, 1])

    c.write(join(test_dir, 'FourPlusSix_conf1.mol'), 0)
    c.write(join(test_dir, 'FourPlusSix_conf2.mol'), 1)
    c.write(join(test_dir, 'FourPlusSix_conf3.mol'), 2)
    c.write(join(test_dir, 'FourPlusSix_conf4.mol'), 3)

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
