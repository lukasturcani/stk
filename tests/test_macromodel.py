"""
Tests functions which use MacroModel.

These tests are only run when the --macromodel_path pytest option is
used. MacroModel is 3rd party software, it does not come with ``stk``.

"""

import pytest
import sys
import os
from os.path import join
import numpy as np
import stk

macromodel = pytest.mark.skipif(
    all('macromodel' not in x for x in sys.argv),
    reason="Only run when explicitly asked.")


outdir = 'macromodel_tests_output'
if not os.path.exists(outdir):
    os.mkdir(outdir)


@macromodel
def test_restricted_force_field(tmp_cc3, macromodel_path):
    tmp_cc3.write(join(outdir, 'rmm_ff_before.mol'), conformer=0)

    mm = stk.MacroModelForceField(macromodel_path=macromodel_path,
                                  output_dir='rmm_ff',
                                  restricted=True,
                                  minimum_gradient=1)
    mm.optimize(tmp_cc3, conformer=0)
    tmp_cc3.write(join(outdir, 'rmm_ff_after.mol'), conformer=0)


@macromodel
def test_unrestricted_force_field(tmp_cc3, macromodel_path):
    tmp_cc3.write(join(outdir, 'umm_ff_before.mol'),
                  conformer=0)

    mm = stk.MacroModelForceField(macromodel_path=macromodel_path,
                                  output_dir='umm_ff',
                                  restricted=False,
                                  minimum_gradient=1)
    mm.optimize(tmp_cc3, conformer=0)
    tmp_cc3.write(join(outdir, 'umm_ff_after.mol'), conformer=0)


@macromodel
def test_both_force_field(tmp_cc3, macromodel_path):
    tmp_cc3.write(join(outdir, 'bmm_ff_before.mol'),
                  conformer=0)

    mm = stk.MacroModelForceField(macromodel_path=macromodel_path,
                                  output_dir='bmm_ff',
                                  restricted='both',
                                  minimum_gradient=1)
    mm.optimize(tmp_cc3, conformer=0)
    tmp_cc3.write(join(outdir, 'bmm_ff_after.mol'), conformer=0)


@macromodel
def test_restricted_md(tmp_cc3, macromodel_path):
    tmp_cc3.write(join(outdir, 'rmm_md_before.mol'), conformer=0)

    # Freeze one of the bonders.
    bonder = tmp_cc3.func_groups[0].bonder_ids[0]
    restricted_bonds = []
    for neighbor in tmp_cc3.mol.GetAtomWithIdx(bonder).GetNeighbors():
        restricted_bonds.append(frozenset((bonder, neighbor.GetIdx())))

    mm = stk.MacroModelMD(macromodel_path=macromodel_path,
                          output_dir='rmm_md',
                          minimum_gradient=1,
                          simulation_time=20,
                          eq_time=2,
                          conformers=2,
                          restricted_bonds=restricted_bonds)
    mm.optimize(tmp_cc3, conformer=0)
    tmp_cc3.write(join(outdir, 'rmm_md_after.mol'), conformer=0)


@macromodel
def test_unrestricted_md(tmp_cc3, macromodel_path):
    tmp_cc3.write(join(outdir, 'umm_md_before.mol'), conformer=0)

    mm = stk.MacroModelMD(macromodel_path=macromodel_path,
                          output_dir='umm_md',
                          minimum_gradient=1,
                          simulation_time=20,
                          eq_time=2,
                          conformers=2)
    mm.optimize(tmp_cc3, conformer=0)
    tmp_cc3.write(join(outdir, 'umm_md_after.mol'), conformer=0)


@macromodel
def test_energy(amine2, macromodel_path):
    mm = stk.MacroModelEnergy(macromodel_path, 'energy_calc')
    a = mm.energy(amine2)
    assert np.allclose(a=a,
                       b=49.0655,
                       atol=1e-2)
