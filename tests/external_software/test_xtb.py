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
    all('xtb_path' not in x for x in sys.argv),
    reason="Only run when explicitly asked.")


@xtb
def test_xtb_negfreq(tmp_polymer, xtb_path):
    # XTB requires an embedding before working.
    etkdg = stk.ETKDG()
    etkdg.optimize(tmp_polymer)

    out_dir = 'gfnxtb_NF_energy'
    energy = stk.XTBEnergy(
        xtb_path=xtb_path,
        output_dir=join(odir, out_dir),
        unlimited_memory=False
    )

    # Calculate the initial energy.
    init_energy = energy.get_energy(tmp_polymer)

    # Run a low criteria optimization.
    out_dir = 'gfnxtb_NF_crude_opt'
    tmp_polymer.write(join(odir, 'gfnxtb_opt_before.mol'))
    opt_lowtol = stk.XTB(
        xtb_path=xtb_path,
        output_dir=join(odir, out_dir),
        unlimited_memory=False,
        opt_level='crude',
        max_runs=2,
        calculate_hessian=True
    )
    opt_lowtol.optimize(tmp_polymer)
    tmp_polymer.write(join(odir, 'gfnxtb_opt_after.mol'))
    # Check that the low criteria calculation is incomplete.
    assert tmp_polymer in opt_lowtol.incomplete
    # Check for restart file.
    assert os.path.isfile(join(join(odir, out_dir), 'xtbhess.coord'))
    # Check number of output files is more than one because the
    # optimisation is crude and max_runs = 2.
    assert len(glob.glob(f'{join(odir, out_dir)}/*.output')) > 1
    # Optimized structure has lower energy than initial structure.
    new_energy = energy.get_energy(tmp_polymer)
    assert new_energy < init_energy

    # Run high criteria optimization.
    out_dir = 'gfnxtb_NF_extreme_opt'
    tmp_polymer.write(join(odir, 'gfnxtb_opt_before.mol'))
    opt_hightol = stk.XTB(
        xtb_path=xtb_path,
        output_dir=join(odir, out_dir),
        unlimited_memory=False,
        opt_level='extreme',
        max_runs=1,
        calculate_hessian=True
    )
    opt_hightol.optimize(tmp_polymer)
    tmp_polymer.write(
        join(join(odir, out_dir), 'gfnxtb_opt_after.mol')
    )
    # Check that the calculation is complete.
    assert tmp_polymer not in opt_hightol.incomplete
    # Check that there is no restart file.
    assert (
        not os.path.isfile(join(join(odir, out_dir), 'xtbhess.coord'))
    )

    # Check that the calculation was run only once.
    assert len(glob.glob(f'{join(odir, out_dir)}/*.output')) == 1

    # Optimized structure has lower energy than initial structure.
    new_energy = energy.get_energy(tmp_polymer)
    assert new_energy < init_energy


@xtb
def test_xtb_solvent_charge_uhf(tmp_polymer, xtb_path):
    init_dir = os.getcwd()
    # XTB requires an embedding before working.
    etkdg = stk.ETKDG()
    etkdg.optimize(tmp_polymer)

    # Test that the energies calculated using xTB with implicit
    # solvation, non zero charge, or non zero num_unpaired_electrons is
    # different to the default case.
    out_dir = 'gfnxtb_energy'
    energy = stk.XTBEnergy(
        xtb_path=xtb_path,
        output_dir=join(odir, out_dir),
        unlimited_memory=True
    )
    init_energy = energy.get_energy(tmp_polymer)

    # Check that directory movement has worked.
    assert os.getcwd() == init_dir

    out_dir = 'gfnxtb_solv_energy'
    solvent = stk.XTBEnergy(
        xtb_path=xtb_path,
        output_dir=join(odir, out_dir),
        unlimited_memory=True,
        solvent='h2o'
    )
    solv_energy = solvent.get_energy(tmp_polymer)
    assert solv_energy != init_energy

    # Check that directory movement has worked.
    assert os.getcwd() == init_dir

    out_dir = 'gfnxtb_charge_energy'
    charge = stk.XTBEnergy(
        xtb_path=xtb_path,
        output_dir=join(odir, out_dir),
        unlimited_memory=True,
        charge=-1
    )
    charge_energy = charge.get_energy(tmp_polymer)
    assert charge_energy != init_energy

    # Check that directory movement has worked.
    assert os.getcwd() == init_dir

    out_dir = 'gfnxtb_multi_energy'
    multi = stk.XTBEnergy(
        xtb_path=xtb_path,
        output_dir=join(odir, out_dir),
        unlimited_memory=True,
        num_unpaired_electrons=2
    )
    multi_energy = multi.get_energy(tmp_polymer)
    assert multi_energy != init_energy

    # Check that directory movement has worked.
    assert os.getcwd() == init_dir


@xtb
def test_xtb_properties(tmp_polymer, xtb_path):
    # XTB requires an embedding before working.
    etkdg = stk.ETKDG()
    etkdg.optimize(tmp_polymer)

    # Calculation of the Hessian requires a well optimized structure.
    gfnxtb = stk.XTB(
        xtb_path,
        output_dir=join(odir, 'gfnxtb_opt'),
        unlimited_memory=False,
        opt_level='extreme',
        max_runs=1,
        calculate_hessian=False,
        num_cores=2,
    )
    gfnxtb.optimize(tmp_polymer)
    xtb = stk.XTBEnergy(
        xtb_path=xtb_path,
        output_dir=join(odir, 'gfnxtb_ey'),
        unlimited_memory=False,
        calculate_free_energy=True
    )
    total_energy = xtb.get_energy(tmp_polymer)

    total_free_energy = xtb.total_free_energies[tmp_polymer]
    assert total_energy != total_free_energy
    assert isinstance(xtb.total_energies[tmp_polymer], float)
    assert isinstance(xtb.total_free_energies[tmp_polymer], float)
    assert isinstance(xtb.frequencies[tmp_polymer], list)
    assert len(xtb.frequencies[tmp_polymer]) > 0
    for i in xtb.frequencies[tmp_polymer]:
        assert isinstance(i, float)
    assert isinstance(xtb.homo_lumo_gaps[tmp_polymer], float)
    assert isinstance(xtb.fermi_levels[tmp_polymer], float)

    assert isinstance(xtb.qonly_dipole_moments[tmp_polymer], list)
    assert len(xtb.qonly_dipole_moments[tmp_polymer]) == 3
    for i in xtb.qonly_dipole_moments[tmp_polymer]:
        assert isinstance(i, float)

    assert isinstance(xtb.full_dipole_moments[tmp_polymer], list)
    assert len(xtb.full_dipole_moments[tmp_polymer]) == 4
    for i in xtb.full_dipole_moments[tmp_polymer]:
        assert isinstance(i, float)

    assert isinstance(xtb.qonly_quadrupole_moments[tmp_polymer], list)
    assert len(xtb.qonly_quadrupole_moments[tmp_polymer]) == 6
    for i in xtb.qonly_quadrupole_moments[tmp_polymer]:
        assert isinstance(i, float)

    assert isinstance(xtb.qdip_quadrupole_moments[tmp_polymer], list)
    assert len(xtb.qdip_quadrupole_moments[tmp_polymer]) == 6
    for i in xtb.qdip_quadrupole_moments[tmp_polymer]:
        assert isinstance(i, float)
    assert isinstance(xtb.full_quadrupole_moments[tmp_polymer], list)
    assert len(xtb.full_quadrupole_moments[tmp_polymer]) == 6
    for i in xtb.full_quadrupole_moments[tmp_polymer]:
        assert isinstance(i, float)

    assert isinstance(xtb.homo_lumo_orbitals[tmp_polymer], dict)
    assert (
        isinstance(xtb.homo_lumo_orbitals[tmp_polymer]['HOMO'], list)
    )
    assert len(xtb.homo_lumo_orbitals[tmp_polymer]['HOMO']) == 3
    for i in xtb.homo_lumo_orbitals[tmp_polymer]['HOMO']:
        assert isinstance(i, float) or isinstance(i, int)

    assert isinstance(xtb.homo_lumo_orbitals[tmp_polymer], dict)
    assert (
        isinstance(xtb.homo_lumo_orbitals[tmp_polymer]['LUMO'], list)
    )
    assert len(xtb.homo_lumo_orbitals[tmp_polymer]['LUMO']) == 3
    for i in xtb.homo_lumo_orbitals[tmp_polymer]['LUMO']:
        assert isinstance(i, float) or isinstance(i, int)


@xtb
def test_valid_solvent():
    gfn_1 = 1
    gfn_2 = 2
    valid = 'acetone'
    valid1 = 'benzene'
    valid2 = 'dmf'
    invalid = 'andrewtarziawrotethis'
    assert stk.is_valid_xtb_solvent(solvent=valid, gfn_version=gfn_1)
    assert stk.is_valid_xtb_solvent(solvent=valid, gfn_version=gfn_2)
    assert stk.is_valid_xtb_solvent(solvent=valid1, gfn_version=gfn_1)
    assert stk.is_valid_xtb_solvent(solvent=valid2, gfn_version=gfn_2)
    assert (
        stk.is_valid_xtb_solvent(solvent=valid1, gfn_version=gfn_2)
        is False
    )
    assert (
        stk.is_valid_xtb_solvent(solvent=valid2, gfn_version=gfn_1)
        is False
    )
    assert (
        stk.is_valid_xtb_solvent(solvent=invalid, gfn_version=gfn_1)
        is False
    )
    assert (
        stk.is_valid_xtb_solvent(solvent=invalid, gfn_version=gfn_2)
        is False
    )
