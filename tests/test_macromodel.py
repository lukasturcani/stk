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
import rdkit.Chem.AllChem as rdkit

macromodel = pytest.mark.skipif(
    all('macromodel' not in x for x in sys.argv),
    reason="Only run when explicitly asked.")


outdir = 'macromodel_tests_output'
if not os.path.exists(outdir):
    os.mkdir(outdir)


@macromodel
def test_restricted_force_field(tmp_cc3, macromodel_path):
    bonder_ids = {
        bid for fg in tmp_cc3.func_groups for bid in fg.bonder_ids
    }

    # Save all bond lengths, angles and torsions which are restricted.
    dims_before = {}
    for bond in tmp_cc3.mol.GetBonds():
        a1, a2 = atoms = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        if a1 not in bonder_ids and a2 not in bonder_ids:
            bond_length = tmp_cc3.atom_distance(a1, a2, conformer=0)
            dims_before[(min(atoms), max(atoms))] = bond_length

    conf = tmp_cc3.mol.GetConformer(0)
    bond_angles = rdkit.FindAllPathsOfLengthN(
        mol=tmp_cc3.mol,
        length=3,
        useBonds=False,
        useHs=True
    )
    for atoms in bond_angles:
        dims_before[tuple(atoms)] = rdkit.GetAngleDeg(conf, *atoms)

    torsion_angles = rdkit.FindAllPathsOfLengthN(
        mol=tmp_cc3.mol,
        length=4,
        useBonds=False,
        useHs=True
    )
    for atoms in torsion_angles:
        dims_before[tuple(atoms)] = rdkit.GetDihedralDeg(conf, *atoms)

    tmp_cc3.write(join(outdir, 'rmm_ff_before.mol'), conformer=0)

    mm = stk.MacroModelForceField(
        macromodel_path=macromodel_path,
        output_dir='rmm_ff',
        restricted=True,
        minimum_gradient=1
    )
    mm.optimize(tmp_cc3, conformer=0)
    tmp_cc3.write(join(outdir, 'rmm_ff_after.mol'), conformer=0)

    # Save all bond lengths, angles and torsions which are restricted.
    dims_after = {}
    for bond in tmp_cc3.mol.GetBonds():
        a1, a2 = atoms = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        if a1 not in bonder_ids and a2 not in bonder_ids:
            bond_length = tmp_cc3.atom_distance(a1, a2, conformer=0)
            dims_after[(min(atoms), max(atoms))] = bond_length

    conf = tmp_cc3.mol.GetConformer(0)
    bond_angles = rdkit.FindAllPathsOfLengthN(
        mol=tmp_cc3.mol,
        length=3,
        useBonds=False,
        useHs=True
    )
    for atoms in bond_angles:
        dims_after[tuple(atoms)] = rdkit.GetAngleDeg(conf, *atoms)

    torsion_angles = rdkit.FindAllPathsOfLengthN(
        mol=tmp_cc3.mol,
        length=4,
        useBonds=False,
        useHs=True
    )
    for atoms in torsion_angles:
        dims_after[tuple(atoms)] = rdkit.GetDihedralDeg(conf, *atoms)

    assert dims_before == dims_after


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
def test_restricted_md(tmp_cc3, macromodel_path):
    tmp_cc3.write(join(outdir, 'rmm_md_before.mol'), conformer=0)

    # Freeze all bond lengths for one of the bonders.
    bonder = tmp_cc3.func_groups[0].bonder_ids[0]
    restricted_bonds = []
    for neighbor in tmp_cc3.mol.GetAtomWithIdx(bonder).GetNeighbors():
        restricted_bonds.append(frozenset((bonder, neighbor.GetIdx())))

    # Freeze some bond angle.
    bond_angles = rdkit.FindAllPathsOfLengthN(
        mol=tmp_cc3.mol,
        length=3,
        useBonds=False,
        useHs=True
    )
    restricted_bond_angles = [frozenset(bond_angles[0])]

    # Freeze some torsional angle.
    torsion_angles = rdkit.FindAllPathsOfLengthN(
        mol=tmp_cc3.mol,
        length=4,
        useBonds=False,
        useHs=True
    )
    restricted_torsional_angles = [frozenset(torsion_angles[0])]

    # Save all bond lengths, angles and torsions which are restricted.
    dims_before = {}
    for atoms in restricted_bonds:
        bond_length = tmp_cc3.atom_distance(*atoms, conformer=0)
        dims_before[atoms] = bond_length

    conf = tmp_cc3.mol.GetConformer(0)
    for atoms in restricted_bond_angles:
        dims_before[atoms] = rdkit.GetAngleDeg(conf, *atoms)

    for atoms in restricted_torsional_angles:
        dims_before[atoms] = rdkit.GetDihedralDeg(conf, *atoms)

    mm = stk.MacroModelMD(
        macromodel_path=macromodel_path,
        output_dir='rmm_md',
        minimum_gradient=1,
        simulation_time=20,
        eq_time=2,
        conformers=2,
        time_step=0.1,
        restricted_bonds=restricted_bonds,
        restricted_bond_angles=restricted_bond_angles,
        restricted_torsional_angles=restricted_torsional_angles
    )
    mm.optimize(tmp_cc3, conformer=0)
    tmp_cc3.write(join(outdir, 'rmm_md_after.mol'), conformer=0)

    # Save all bond lengths, angles and torsions which are restricted.
    dims_after = {}
    for atoms in restricted_bonds:
        bond_length = tmp_cc3.atom_distance(*atoms, conformer=0)
        dims_after[atoms] = bond_length

    conf = tmp_cc3.mol.GetConformer(0)
    for atoms in restricted_bond_angles:
        dims_after[atoms] = rdkit.GetAngleDeg(conf, *atoms)

    for atoms in restricted_torsional_angles:
        dims_after[atoms] = rdkit.GetDihedralDeg(conf, *atoms)

    assert dims_before == dims_after


@macromodel
def test_unrestricted_md(tmp_cc3, macromodel_path):
    tmp_cc3.write(join(outdir, 'umm_md_before.mol'), conformer=0)

    mm = stk.MacroModelMD(macromodel_path=macromodel_path,
                          output_dir='umm_md',
                          minimum_gradient=1,
                          simulation_time=20,
                          eq_time=2,
                          conformers=2,
                          time_step=0.1)
    mm.optimize(tmp_cc3, conformer=0)
    tmp_cc3.write(join(outdir, 'umm_md_after.mol'), conformer=0)


@macromodel
def test_energy(amine2, macromodel_path):
    mm = stk.MacroModelEnergy(macromodel_path, 'energy_calc')
    a = mm.energy(amine2)
    assert np.allclose(a=a,
                       b=49.0655,
                       atol=1e-2)
