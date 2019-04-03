import os
from os.path import join
import stk


if not os.path.exists('fg_tests_output'):
    os.mkdir('fg_tests_output')


def test_bi_fg_bb(diol2, difluorene_dibromine):
    p = stk.Polymer([diol2, difluorene_dibromine],
                    stk.Linear('ABAB', [1, 1, 0, 0], 4))
    p.write(join('fg_tests_output/diol_difluorene_dibromine.mol'))


def test_diol_with_difluorene(diol2, difluorene2):
    p = stk.Polymer([diol2, difluorene2],
                    stk.Linear('ABAB', [1, 1, 0, 0], 4))
    p.write(join('fg_tests_output/diol_difluorene.mol'))


def test_boronic_acid_with_diol(boronic_acid2, diol2):
    p = stk.Polymer([boronic_acid2, diol2],
                    stk.Linear('ABAB', [1, 1, 0, 0], 4))
    p.write(join('fg_tests_output/boronic_acid_diol.mol'))


def test_phenyl_with_phenyl_amine(phenyl_amine):
    p = stk.Polymer([phenyl_amine], stk.Linear('A', [1], 8))
    p.write(join('fg_tests_output/phenyl_amine_phenyl_amine.mol'))
