import os
from os.path import join
import stk

test_dir = 'functional_group_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


def test_remove_deleters(tmp_fg):
    tmp_fg.remove_deleters([1, 2, 15, 55, 100, 300])
    assert tmp_fg.atom_ids == (8, 1, 2, 40, 3, 29)
    assert tmp_fg.bonder_ids == (1, 29, 8)
    assert tmp_fg.deleter_ids == (3, )


def test_shifted_fg(fg):
    shifted = fg.shifted_fg(fg.id+1, 20)
    assert shifted is not fg
    assert shifted.id == fg.id+1

    for a1, a2 in zip(fg.atom_ids, shifted.atom_ids):
        assert a1 + 20 == a2

    for a1, a2 in zip(fg.bonder_ids, shifted.bonder_ids):
        assert a1 + 20 == a2

    for a1, a2 in zip(fg.bonder_ids, shifted.bonder_ids):
        assert a1 + 20 == a2


def test_bi_fg_bb(diol2, difluorene_dibromine):
    p = stk.Polymer([diol2, difluorene_dibromine],
                    stk.Linear('ABAB', [1, 1, 0, 0], 4))
    p.write(join(test_dir, 'diol_difluorene_dibromine.mol'))


def test_diol_with_difluorene(diol2, difluorene2):
    p = stk.Polymer([diol2, difluorene2],
                    stk.Linear('ABAB', [1, 1, 0, 0], 4))
    p.write(join(test_dir, 'diol_difluorene.mol'))


def test_boronic_acid_with_diol(boronic_acid2, diol2):
    p = stk.Polymer([boronic_acid2, diol2],
                    stk.Linear('ABAB', [1, 1, 0, 0], 4))
    p.write(join(test_dir, 'boronic_acid_diol.mol'))


def test_phenyl_with_phenyl_amine(phenyl_amine):
    p = stk.Polymer([phenyl_amine], stk.Linear('A', [1], 8))
    p.write(join(test_dir, 'phenyl_amine_phenyl_amine.mol'))
