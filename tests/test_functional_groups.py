import os
from os.path import join
import stk

test_dir = 'functional_group_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


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


def test_phenyl_with_ring_amine(ring_amine):
    p = stk.Polymer([ring_amine], stk.Linear('A', [1], 8))
    p.write(join(test_dir, 'ring_amine_ring_amine.mol'))
