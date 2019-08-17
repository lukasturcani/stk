"""
Tests functions which use MacroModel.

These tests are only run when the --macromodel_path pytest option is
used. MacroModel is 3rd party software, it does not come with ``stk``.

"""

import pytest
import sys
import os
from os.path import join
import stk

macromodel = pytest.mark.skipif(
    all('macromodel' not in x for x in sys.argv),
    reason="Only run when explicitly asked.")


test_dir = 'macromodel_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


@macromodel
def test_restricted_force_field(tmp_tetrahedron, macromodel_path):

    mm_energy = stk.MacroModelEnergy(macromodel_path, force_field=16)
    init_energy = mm_energy.get_energy(tmp_tetrahedron)
    tmp_tetrahedron.write(join(test_dir, 'rmm_ff_before.mol'))

    mm = stk.MacroModelForceField(
        macromodel_path=macromodel_path,
        output_dir='rmm_ff',
        restricted=True,
        minimum_gradient=1,
        force_field=16
    )
    mm.optimize(tmp_tetrahedron)
    tmp_tetrahedron.write(join(test_dir, 'rmm_ff_after.mol'))

    assert mm_energy.get_energy(tmp_tetrahedron) < init_energy


@macromodel
def test_unrestricted_force_field(tmp_tetrahedron, macromodel_path):

    mm_energy = stk.MacroModelEnergy(macromodel_path, force_field=16)
    init_energy = mm_energy.get_energy(tmp_tetrahedron)
    tmp_tetrahedron.write(join(test_dir, 'umm_ff_before.mol'))

    mm = stk.MacroModelForceField(
        macromodel_path=macromodel_path,
        output_dir='umm_ff',
        restricted=False,
        minimum_gradient=1,
        force_field=16
    )
    mm.optimize(tmp_tetrahedron)
    tmp_tetrahedron.write(join(test_dir, 'umm_ff_after.mol'))

    assert mm_energy.get_energy(tmp_tetrahedron) < init_energy


@macromodel
def test_restricted_md(tmp_tetrahedron, macromodel_path):

    mm_energy = stk.MacroModelEnergy(macromodel_path, force_field=16)
    init_energy = mm_energy.get_energy(tmp_tetrahedron)
    tmp_tetrahedron.write(join(test_dir, 'rmm_md_before.mol'))

    # Freeze one of the bonders.
    bonder = tmp_tetrahedron.func_groups[0].bonders[0].id
    rdkit_mol = tmp_tetrahedron.to_rdkit_mol()
    restricted_bonds = []
    for neighbor in rdkit_mol.GetAtomWithIdx(bonder).GetNeighbors():
        restricted_bonds.append(frozenset((bonder, neighbor.GetIdx())))

    mm = stk.MacroModelMD(
        macromodel_path=macromodel_path,
        output_dir='rmm_md',
        minimum_gradient=1,
        simulation_time=20,
        eq_time=2,
        conformers=2,
        time_step=0.1,
        force_field=16,
        restricted_bonds=restricted_bonds
    )
    mm.optimize(tmp_tetrahedron)
    tmp_tetrahedron.write(join(test_dir, 'rmm_md_after.mol'))

    assert mm_energy.get_energy(tmp_tetrahedron) < init_energy


@macromodel
def test_unrestricted_md(tmp_tetrahedron, macromodel_path):

    mm_energy = stk.MacroModelEnergy(macromodel_path, force_field=16)
    init_energy = mm_energy.get_energy(tmp_tetrahedron)
    tmp_tetrahedron.write(join(test_dir, 'umm_md_before.mol'))

    mm = stk.MacroModelMD(
        macromodel_path=macromodel_path,
        output_dir='umm_md',
        minimum_gradient=1,
        simulation_time=20,
        eq_time=2,
        conformers=2,
        time_step=0.1,
        force_field=16
    )
    mm.optimize(tmp_tetrahedron)
    tmp_tetrahedron.write(join(test_dir, 'umm_md_after.mol'))

    assert mm_energy.get_energy(tmp_tetrahedron) < init_energy


@macromodel
def test_energy(amine2, amine2_conf1, macromodel_path):
    mm = stk.MacroModelEnergy(macromodel_path, 'energy_calc')
    assert mm.get_energy(amine2) < mm.get_energy(amine2_conf1)


def test_forcefield_com_exceptions():
    with pytest.raises(stk.MacroModelInputError):
        stk.MacroModelForceField(
            macromodel_path='dummy_path',
            maximum_iterations=1000000
        )

    with pytest.raises(stk.MacroModelInputError):
        stk.MacroModelForceField(
            macromodel_path='dummy_path',
            minimum_gradient=0.00001
        )


def test_md_com_exceptions(amine2):
    with pytest.raises(stk.MacroModelInputError):
        mm = stk.MacroModelMD(
            macromodel_path='dummy_path',
            conformers=10000
        )

    with pytest.raises(stk.MacroModelInputError):
        mm = stk.MacroModelMD(
            macromodel_path='dummy_path',
            simulation_time=1000000
        )

    with pytest.raises(stk.MacroModelInputError):
        mm = stk.MacroModelMD(
            macromodel_path='dummy_path',
            time_step=100000
        )

    with pytest.raises(stk.MacroModelInputError):
        mm = stk.MacroModelMD(
            macromodel_path='dummy_path',
            eq_time=1000000
        )

    with pytest.raises(stk.MacroModelInputError):
        mm = stk.MacroModelMD(
            macromodel_path='dummy_path',
            maximum_iterations=1000000
        )

    with pytest.raises(stk.MacroModelInputError):
        mm = stk.MacroModelMD(
            macromodel_path='dummy_path',
            minimum_gradient=0.00001
        )

    mm = stk.MacroModelMD(
        macromodel_path='dummy_path',
        simulation_time=100000,
        eq_time=100000
    )

    mm._generate_com(amine2, join(test_dir, 'com_test'))
    with open(join(test_dir, 'com_test.com'), 'r') as o:
        comfile = o.read().splitlines()
        expect1 = (
            ' MDYN       0      0      0      0     1.0000 -10'
            '00.0000   300.0000     0.0000'
        )
        expect2 = (
            ' MDYN       1      0      0      0     1.0000 -10'
            '00.0000   300.0000     0.0000'
        )
        assert comfile[5] == expect1
        assert comfile[7] == expect2
