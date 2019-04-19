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
    always_raiser = stk.RaisingOptimizer(fn=etkdg.optimize,
                                         fail_chance=1)
    with pytest.raises(stk.RaisingOptimizerError):
        always_raiser.optimize(tmp_polymer)

    never_raiser = stk.RaisingOptimizer(fn=etkdg.optimize,
                                        fail_chance=0)
    tmp_polymer.write(join(odir, 'raising_optimizer_before.mol'))
    never_raiser.optimize(tmp_polymer)
    tmp_polymer.write(join(odir, 'raising_optimizer_after.mol'))


def test_rdkit_force_field(tmp_polymer):
    tmp_polymer.write(join(odir, 'rdkit_force_field_before.mol'))
    mmff = stk.RDKitForceField(rdkit.MMFFOptimizeMolecule)
    mmff.optimize(tmp_polymer)
    tmp_polymer.write(join(odir, 'rdkit_force_field_after.mol'))


def test_rdkit_embedder(tmp_polymer):
    tmp_polymer.write(join(odir, 'rdkit_embed_before.mol'))
    etkdg = stk.RDKitEmbedder(rdkit.ETKDG())
    etkdg.optimize(tmp_polymer)
    tmp_polymer.write(join(odir, 'rdkit_embed_after.mol'))


def test_optimizer_pipeline(tmp_polymer):
    tmp_polymer.write(join(odir, 'pipeline_before.mol'))
    etkdg = stk.RDKitEmbedder(rdkit.ETKDG())
    mmff = stk.RDKitForceField(rdkit.MMFFOptimizeMolecule)
    pipeline = stk.OptimizerPipeline(etkdg, mmff)
    pipeline.optimize(tmp_polymer)
    tmp_polymer.write(join(odir, 'pipeline_after.mol'))


def test_optimizer_skipping(tmp_polymer):
    etkdg = stk.RDKitEmbedder(rdkit.ETKDG())
    skipper = stk.RDKitEmbedder(rdkit.ETKDG(), skip_optimized=True)

    etkdg.optimize(tmp_polymer)
    etkdg.optimize(tmp_polymer)
    skipper.optimize(tmp_polymer)


def test_cage_pipeline(tmp_cc3, tmp_cage):
    mmff = stk.RDKitForceField(rdkit.MMFFOptimizeMolecule)
    etkdg = stk.RDKitEmbedder(rdkit.ETKDG())
    pipeline = stk.CageOptimizerPipeline(etkdg, mmff)
    pipeline.optimize(tmp_cage)
    pipeline.optimize(tmp_cc3)
