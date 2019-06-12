import stk
from os.path import join
import os
import pytest
import glob
import sys

odir = 'optimizer_tests_output'
if not os.path.exists(odir):
    os.mkdir(odir)

xtb = pytest.mark.skipif(
    all('xtb' not in x for x in sys.argv),
    reason="Only run when explicitly asked.")


@xtb
def test_xtb_negfreq(tmp_polymer, xtb_path):
    init_dir = os.getcwd()
    # XTB  requires an embedding before working.
    etkdg = stk.ETKDG()
    etkdg.optimize(tmp_polymer)

    conformer = tmp_polymer.mol.GetConformer(-1).GetId()
    id = (tmp_polymer, conformer)

    out_dir = 'gfnxtb_NF_energy'
    energy = stk.XTBEnergy(
        xtb_path=xtb_path,
        output_dir=join(odir, out_dir),
        mem_ulimit=True
    )

    # initial energy
    init_energy = energy.energy(tmp_polymer)

    # run low criteria optimization that will loop until negative frequencies
    # are removed
    out_dir = 'gfnxtb_NF_crude_opt'
    tmp_polymer.write(join(odir, 'gfnxtb_opt_before.mol'))
    opt_lowtol = stk.XTB(
        xtb_path=xtb_path,
        output_dir=join(odir, out_dir),
        mem_ulimit=True,
        opt_level='crude',
        max_runs=2
    )
    opt_lowtol.optimize(tmp_polymer)
    tmp_polymer.write(join(odir, 'gfnxtb_opt_after.mol'))
    # check incomplete flag
    assert id in opt_lowtol.incomplete
    # check for restart file
    assert os.path.isfile(join(join(odir, out_dir), 'xtbhess.coord'))
    # check number of output files is more than one
    if len(glob.glob(f'*.output')) == 1:
        assert False
    # optimized structure has lower energy than initial structure
    new_energy = energy.energy(tmp_polymer)
    assert new_energy < init_energy

    # high criteria optimization
    out_dir = 'gfnxtb_NF_extreme_opt'
    tmp_polymer.write(join(odir, 'gfnxtb_opt_before.mol'))
    opt_hightol = stk.XTB(
        xtb_path=xtb_path,
        output_dir=join(odir, out_dir),
        mem_ulimit=True,
        opt_level='extreme',
        max_runs=1
    )
    opt_hightol.optimize(tmp_polymer)
    tmp_polymer.write(join(join(odir, out_dir), 'gfnxtb_opt_after.mol'))
    # check incomplete flag
    assert id not in opt_hightol.incomplete
    # check for restart file
    if os.path.isfile(join(join(odir, out_dir), 'xtbhess.coord')):
        assert False
    # check number of output files is more than one
    if len(glob.glob(f'*.output')) == 1:
        assert True
    # optimized structure has lower energy than initial structure
    new_energy = energy.energy(tmp_polymer)
    assert new_energy < init_energy
    # check directory moving worked
    assert os.getcwd() == init_dir


def test_xtb_solvent_charge_unpaired_electrons(tmp_polymer, xtb_path):
    # XTB  requires an embedding before working.
    etkdg = stk.ETKDG()
    etkdg.optimize(tmp_polymer)

    # test that the energies with implicit solvation, charge != 0 and non zero
    # unpaired_electrons
    out_dir = 'gfnxtb_energy'
    energy = stk.XTBEnergy(
        xtb_path=xtb_path,
        output_dir=join(odir, out_dir),
        mem_ulimit=True
    )
    # initial energy
    init_energy = energy.energy(tmp_polymer)

    out_dir = 'gfnxtb_solv_energy'
    solvent = stk.XTBEnergy(
        xtb_path=xtb_path,
        output_dir=join(odir, out_dir),
        mem_ulimit=True,
        solvent='h2o'
    )
    # initial energy
    solv_energy = solvent.energy(tmp_polymer)
    assert solv_energy != init_energy

    out_dir = 'gfnxtb_charge_energy'
    charge = stk.XTBEnergy(
        xtb_path=xtb_path,
        output_dir=join(odir, out_dir),
        mem_ulimit=True,
        charge=-1
    )
    # initial energy
    charge_energy = charge.energy(tmp_polymer)
    assert charge_energy != init_energy

    out_dir = 'gfnxtb_multi_energy'
    multi = stk.XTBEnergy(
        xtb_path=xtb_path,
        output_dir=join(odir, out_dir),
        mem_ulimit=True,
        unpaired_electrons=2
    )
    # initial energy
    multi_energy = multi.energy(tmp_polymer)
    assert multi_energy != init_energy


def test_xtb_properties(tmp_polymer, xtb_path):
    # XTB requires an embedding before working.
    etkdg = stk.ETKDG()
    etkdg.optimize(tmp_polymer)

    conformer = tmp_polymer.mol.GetConformer(-1).GetId()

    # hessian requires well optimized structure
    gfnxtb = stk.XTB(
        xtb_path,
        output_dir=join(odir, 'gfnxtb_opt'),
        mem_ulimit=True,
        opt_level='extreme',
        max_runs=None,
        num_cores=2,
    )
    gfnxtb.optimize(tmp_polymer)
    EC = stk.XTBFreeEnergy(
        xtb_path=xtb_path,
        output_dir=join(odir, 'gfnxtb_ey'),
        mem_ulimit=True
    )
    total_energy = EC.energy(tmp_polymer)
    id = (tmp_polymer, conformer)

    total_free_energy = EC.total_free_energies[id]
    assert total_energy != total_free_energy
    # check directory moving worked
    assert isinstance(EC.total_energies[id], float)
    assert isinstance(EC.total_free_energies[id], float)
    assert isinstance(EC.frequencies[id], list)
    assert len(EC.frequencies[id]) > 0
    for i in EC.frequencies[id]:
        assert isinstance(i, float)
    assert isinstance(EC.homo_lumo_gaps[id], float)
    assert isinstance(EC.fermi_levels[id], float)

    assert isinstance(EC.qonly_dipole_moments[id], list)
    assert len(EC.qonly_dipole_moments[id]) == 3
    for i in EC.qonly_dipole_moments[id]:
        assert isinstance(i, float)

    assert isinstance(EC.full_dipole_moments[id], list)
    assert len(EC.full_dipole_moments[id]) == 4
    for i in EC.full_dipole_moments[id]:
        assert isinstance(i, float)

    assert isinstance(EC.qonly_quadrupole_moments[id], list)
    assert len(EC.qonly_quadrupole_moments[id]) == 6
    for i in EC.qonly_quadrupole_moments[id]:
        assert isinstance(i, float)

    assert isinstance(EC.qdip_quadrupole_moments[id], list)
    assert len(EC.qdip_quadrupole_moments[id]) == 6
    for i in EC.qdip_quadrupole_moments[id]:
        assert isinstance(i, float)
    assert isinstance(EC.full_quadrupole_moments[id], list)
    assert len(EC.full_quadrupole_moments[id]) == 6
    for i in EC.full_quadrupole_moments[id]:
        assert isinstance(i, float)

    assert isinstance(EC.homo_lumo_orbitals[id], dict)
    assert isinstance(EC.homo_lumo_orbitals[id]['HOMO'], list)
    assert len(EC.homo_lumo_orbitals[id]['HOMO']) == 3
    for i in EC.homo_lumo_orbitals[id]['HOMO']:
        assert isinstance(i, float) or isinstance(i, int)

    assert isinstance(EC.homo_lumo_orbitals[id], dict)
    assert isinstance(EC.homo_lumo_orbitals[id]['LUMO'], list)
    assert len(EC.homo_lumo_orbitals[id]['LUMO']) == 3
    for i in EC.homo_lumo_orbitals[id]['LUMO']:
        assert isinstance(i, float) or isinstance(i, int)


def test_valid_solvent():
    gfn_1 = '1'
    gfn_2 = '2'
    valid = 'acetone'  # valid in both versions
    valid1 = 'benzene'  # valid in GFN 1 only
    valid2 = 'dmf'  # valid in GFN 2 only
    invalid = 'andrewtarziawrotethis'
    assert stk.valid_xtb_solvent(solvent=valid, gfn_version=gfn_1)
    assert stk.valid_xtb_solvent(solvent=valid, gfn_version=gfn_2)
    assert stk.valid_xtb_solvent(solvent=valid1, gfn_version=gfn_1)
    assert stk.valid_xtb_solvent(solvent=valid2, gfn_version=gfn_2)
    try:
        stk.valid_xtb_solvent(solvent=valid1, gfn_version=gfn_2)
        assert False
    except stk.XTBInvalidSolventError:
        assert True
    try:
        stk.valid_xtb_solvent(solvent=valid2, gfn_version=gfn_1)
        assert False
    except stk.XTBInvalidSolventError:
        assert True
    try:
        stk.valid_xtb_solvent(solvent=invalid, gfn_version=gfn_1)
        assert False
    except stk.XTBInvalidSolventError:
        assert True
    try:
        stk.valid_xtb_solvent(solvent=invalid, gfn_version=gfn_2)
        assert False
    except stk.XTBInvalidSolventError:
        assert True
