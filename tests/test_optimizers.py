import stk
from os.path import join
import os
import rdkit.Chem.AllChem as rdkit
import pytest
import sys


odir = 'optimizer_tests_output'
if not os.path.exists(odir):
    os.mkdir(odir)


def test_raising_optimizer(tmp_polymer):
    etkdg = stk.RDKitEmbedder(rdkit.ETKDGv2())
    always_raiser = stk.RaisingOptimizer(optimizer=etkdg,
                                         fail_chance=1)
    with pytest.raises(stk.RaisingOptimizerError):
        always_raiser.optimize(tmp_polymer)

    never_raiser = stk.RaisingOptimizer(optimizer=etkdg,
                                        fail_chance=0)
    tmp_polymer.write(join(odir, 'raising_optimizer_before.mol'))
    never_raiser.optimize(tmp_polymer)
    tmp_polymer.write(join(odir, 'raising_optimizer_after.mol'))


def test_mmff(tmp_polymer):
    # If the optimization was successful the energy should be lowered.
    energy_calculator = stk.MMFFEnergy()
    init_energy = energy_calculator.energy(tmp_polymer)

    tmp_polymer.write(join(odir, 'mmff_before.mol'))
    mmff = stk.MMFF()
    mmff.optimize(tmp_polymer)
    tmp_polymer.write(join(odir, 'mmff_after.mol'))

    assert energy_calculator.energy(tmp_polymer) < init_energy


def test_uff(tmp_polymer):
    # If the optimization was successful the energy should be lowered.
    energy_calculator = stk.UFFEnergy()
    init_energy = energy_calculator.energy(tmp_polymer)

    tmp_polymer.write(join(odir, 'uff_before.mol'))
    uff = stk.UFF()
    uff.optimize(tmp_polymer)
    tmp_polymer.write(join(odir, 'uff_after.mol'))

    assert energy_calculator.energy(tmp_polymer) < init_energy


def test_rdkit_embedder(tmp_polymer):
    # If the optimization was successful the energy should be lowered.
    energy_calculator = stk.UFFEnergy()
    init_energy = energy_calculator.energy(tmp_polymer)

    tmp_polymer.write(join(odir, 'rdkit_embed_before.mol'))
    etkdg = stk.RDKitEmbedder(rdkit.ETKDG())
    etkdg.optimize(tmp_polymer)
    tmp_polymer.write(join(odir, 'rdkit_embed_after.mol'))

    assert energy_calculator.energy(tmp_polymer) < init_energy


def test_optimizer_sequence(tmp_polymer):
    # If the optimization was successful the energy should be lowered.
    energy_calculator = stk.MMFFEnergy()
    init_energy = energy_calculator.energy(tmp_polymer)

    tmp_polymer.write(join(odir, 'optimize_sequence_before.mol'))
    etkdg = stk.RDKitEmbedder(rdkit.ETKDG())
    mmff = stk.MMFF()
    sequence = stk.OptimizerSequence(etkdg, mmff)
    sequence.optimize(tmp_polymer)
    tmp_polymer.write(join(odir, 'optimize_sequence_after.mol'))

    assert energy_calculator.energy(tmp_polymer) < init_energy


def test_cache_use(tmp_polymer):
    opt1 = stk.NullOptimizer()
    assert not opt1.cache
    opt1.optimize(tmp_polymer)
    assert not opt1.cache

    opt2 = stk.NullOptimizer(use_cache=True)
    assert not opt2.cache
    opt2.optimize(tmp_polymer)
    assert (tmp_polymer.key, -1) in opt2.cache
    opt2.optimize(tmp_polymer)


def test_cage_optimizer_sequence(tmp_cc3, tmp_cage):
    energy_calculator = stk.MMFFEnergy()
    init_cc3 = energy_calculator.energy(tmp_cc3)
    init_cage = energy_calculator.energy(tmp_cage)

    mmff = stk.MMFF()
    etkdg = stk.RDKitEmbedder(rdkit.ETKDG())

    # CC3 needs an optimizer with mmff only, because using etkdg will
    # increase its energy.
    sequence1 = stk.CageOptimizerSequence(mmff)
    sequence1.optimize(tmp_cc3)
    sequence2 = stk.CageOptimizerSequence(etkdg, mmff)
    sequence2.optimize(tmp_cage)

    # CC3 should have found all windows so energy should be lowered due
    # to optimization.
    assert energy_calculator.energy(tmp_cc3) < init_cc3
    # Cage should have not found all windows so energy should be the
    # same as no optimization happened.
    assert energy_calculator.energy(tmp_cage) == init_cage


def test_try_catch_optimizer(tmp_amine2):
    etkdg = stk.RDKitEmbedder(rdkit.ETKDGv2())
    always_raiser = stk.RaisingOptimizer(optimizer=etkdg,
                                         fail_chance=1)
    success = stk.TryCatchOptimizer(try_optimizer=etkdg,
                                    catch_optimizer=None)
    try_catch = stk.TryCatchOptimizer(try_optimizer=always_raiser,
                                      catch_optimizer=success)
    tmp_amine2.write(
        path=join(odir, 'try_catch_optimizer_before.mol'),
        conformer=1
    )
    try_catch.optimize(tmp_amine2, 1)
    tmp_amine2.write(
        path=join(odir, 'try_catch_optimizer_after.mol'),
        conformer=1
    )


gfnxtb = pytest.mark.skipif(
    all('gfnxtb' not in x for x in sys.argv),
    reason="Only run when explicitly asked.")


@gfnxtb
def test_gfnxtb(tmp_polymer, gfnxtb_path):
    tmp_polymer.write(join(odir, 'gfnxtb_before.mol'))
    gfnxtb = stk.GFNXTB(gfnxtb_path)
    gfnxtb.optimize(tmp_polymer)
    tmp_polymer.write(join(odir, 'gfnxtb_after.mol'))
