import stk
from os.path import join
import os
import pytest


odir = 'optimizer_tests_output'
if not os.path.exists(odir):
    os.mkdir(odir)


def test_raising_optimizer(tmp_polymer):
    mmff = stk.MMFF()
    always_raiser = stk.RaisingOptimizer(
        optimizer=mmff,
        fail_chance=1
    )
    with pytest.raises(stk.RaisingOptimizerError):
        always_raiser.optimize(tmp_polymer)

    never_raiser = stk.RaisingOptimizer(
        optimizer=mmff,
        fail_chance=0
    )
    tmp_polymer.write(join(odir, 'raising_optimizer_before.mol'))

    energy_calculator = stk.MMFFEnergy()
    energy_before = energy_calculator.get_energy(tmp_polymer)

    never_raiser.optimize(tmp_polymer)
    tmp_polymer.write(join(odir, 'raising_optimizer_after.mol'))

    assert energy_before > energy_calculator.get_energy(tmp_polymer)


def test_mmff(tmp_polymer):
    # If the optimization was successful the energy should be lowered.
    energy_calculator = stk.MMFFEnergy()
    init_energy = energy_calculator.get_energy(tmp_polymer)

    tmp_polymer.write(join(odir, 'mmff_before.mol'))
    mmff = stk.MMFF()
    mmff.optimize(tmp_polymer)
    tmp_polymer.write(join(odir, 'mmff_after.mol'))

    assert energy_calculator.get_energy(tmp_polymer) < init_energy


def test_uff(tmp_polymer):
    # If the optimization was successful the energy should be lowered.
    energy_calculator = stk.UFFEnergy()
    init_energy = energy_calculator.get_energy(tmp_polymer)

    tmp_polymer.write(join(odir, 'uff_before.mol'))
    uff = stk.UFF()
    uff.optimize(tmp_polymer)
    tmp_polymer.write(join(odir, 'uff_after.mol'))

    assert energy_calculator.get_energy(tmp_polymer) < init_energy


def test_etkdg(tmp_polymer):
    # If the optimization was successful the energy should be lowered.
    energy_calculator = stk.UFFEnergy()
    init_energy = energy_calculator.get_energy(tmp_polymer)

    tmp_polymer.write(join(odir, 'rdkit_embed_before.mol'))
    etkdg = stk.ETKDG()
    etkdg.optimize(tmp_polymer)
    tmp_polymer.write(join(odir, 'rdkit_embed_after.mol'))

    assert energy_calculator.get_energy(tmp_polymer) < init_energy


def test_optimizer_sequence(tmp_polymer):
    # If the optimization was successful the energy should be lowered.
    energy_calculator = stk.MMFFEnergy()
    init_energy = energy_calculator.get_energy(tmp_polymer)

    tmp_polymer.write(join(odir, 'optimize_sequence_before.mol'))
    etkdg = stk.ETKDG()
    mmff = stk.MMFF()
    sequence = stk.OptimizerSequence(etkdg, mmff)
    sequence.optimize(tmp_polymer)
    tmp_polymer.write(join(odir, 'optimize_sequence_after.mol'))

    assert energy_calculator.get_energy(tmp_polymer) < init_energy


def test_cache_use(tmp_polymer):
    opt1 = stk.NullOptimizer()
    assert not opt1._cache
    opt1.optimize(tmp_polymer)
    assert not opt1._cache

    opt2 = stk.NullOptimizer(use_cache=True)
    assert not opt2._cache
    opt2.optimize(tmp_polymer)
    assert tmp_polymer in opt2._cache


def test_cage_optimizer_sequence(tmp_opt_cc3, tmp_cc3):
    energy_calculator = stk.MMFFEnergy()
    init_opt_cc3 = energy_calculator.get_energy(tmp_opt_cc3)
    init_cc3 = energy_calculator.get_energy(tmp_cc3)

    mmff = stk.MMFF()
    sequence = stk.CageOptimizerSequence(mmff)
    sequence.optimize(tmp_opt_cc3)
    sequence.optimize(tmp_cc3)

    # opt_cc3 should have found all windows so energy should be lowered
    # due to optimization.
    assert energy_calculator.get_energy(tmp_opt_cc3) < init_opt_cc3
    # cc3 should have not found all windows so energy should be the
    # same as if no optimization happened.
    assert energy_calculator.get_energy(tmp_cc3) == init_cc3


def test_try_catch_optimizer(tmp_amine2):
    etkdg = stk.ETKDG()
    always_raiser = stk.RaisingOptimizer(
        optimizer=etkdg,
        fail_chance=1
    )
    success = stk.TryCatchOptimizer(
        try_optimizer=etkdg,
        catch_optimizer=None
    )
    try_catch = stk.TryCatchOptimizer(
        try_optimizer=always_raiser,
        catch_optimizer=success
    )
    tmp_amine2.write(join(odir, 'try_catch_optimizer_before.mol'))
    try_catch.optimize(tmp_amine2)
    tmp_amine2.write(join(odir, 'try_catch_optimizer_after.mol'))
