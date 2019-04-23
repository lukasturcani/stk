import stk
from os.path import join
import os
import rdkit.Chem.AllChem as rdkit
import pytest


odir = 'optimizer_tests_output'
if not os.path.exists(odir):
    os.mkdir(odir)


def test_raising_optimizer(tmp_polymer):
    etkdg = stk.RDKitEmbedder(rdkit.ETKDGv2())
    always_raiser = stk.RaisingOptimizer(optimizer=etkdg.optimize,
                                         fail_chance=1)
    with pytest.raises(stk.RaisingOptimizerError):
        always_raiser.optimize(tmp_polymer)

    never_raiser = stk.RaisingOptimizer(optimizer=etkdg.optimize,
                                        fail_chance=0)
    tmp_polymer.write(join(odir, 'raising_optimizer_before.mol'))
    never_raiser.optimize(tmp_polymer)
    tmp_polymer.write(join(odir, 'raising_optimizer_after.mol'))


def test_mmff(tmp_polymer):
    tmp_polymer.write(join(odir, 'mmff_before.mol'))
    mmff = stk.MMFF()
    mmff.optimize(tmp_polymer)
    tmp_polymer.write(join(odir, 'mmff_after.mol'))


def test_uff(tmp_polymer):
    tmp_polymer.write(join(odir, 'uff_before.mol'))
    mmff = stk.UFF()
    mmff.optimize(tmp_polymer)
    tmp_polymer.write(join(odir, 'uff_after.mol'))


def test_rdkit_embedder(tmp_polymer):
    tmp_polymer.write(join(odir, 'rdkit_embed_before.mol'))
    etkdg = stk.RDKitEmbedder(rdkit.ETKDG())
    etkdg.optimize(tmp_polymer)
    tmp_polymer.write(join(odir, 'rdkit_embed_after.mol'))


def test_optimize_sequence(tmp_polymer):
    tmp_polymer.write(join(odir, 'optimize_sequence_before.mol'))
    etkdg = stk.RDKitEmbedder(rdkit.ETKDG())
    mmff = stk.MMFF()
    pipeline = stk.Sequence(etkdg, mmff)
    pipeline.optimize(tmp_polymer)
    tmp_polymer.write(join(odir, 'optimize_sequence_after.mol'))


def test_cache_use(tmp_polymer):
    etkdg = stk.RDKitEmbedder(rdkit.ETKDG())
    skipper = stk.RDKitEmbedder(rdkit.ETKDG(), skip_optimized=True)

    etkdg.optimize(tmp_polymer)
    etkdg.optimize(tmp_polymer)
    skipper.optimize(tmp_polymer)


def test_cage_optimizer_sequence(tmp_cc3, tmp_cage):
    mmff = stk.MMFF()
    etkdg = stk.RDKitEmbedder(rdkit.ETKDG())
    pipeline = stk.CageOptimizerPipeline(etkdg, mmff)
    pipeline.optimize(tmp_cage)
    pipeline.optimize(tmp_cc3)


def test_try_catch_optimizer(tmp_amine2):
    etkdg = stk.RDKitEmbedder(rdkit.ETKDGv2())
    always_raiser = stk.RaisingOptimizer(optimizer=etkdg,
                                         fail_chance=1)
    success = stk.TryCatchOptimizer(try_optimizer=etkdg,
                                    catch_optimizer=None)
    try_catch = stk.TryCatchOptimizer(try_optimizer=always_raiser,
                                      catch_optimizer=success)
    tmp_amine2.write('try_catch_optimizer_before.mol', 1)
    try_catch.optimize(tmp_amine2, 1)
    tmp_amine2.write('try_catch_optimizer_after.mol', 1)
