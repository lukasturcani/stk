import stk
from os.path import join
import os
import pytest
import sys

odir = 'optimizer_tests_output'
if not os.path.exists(odir):
    os.mkdir(odir)


def amine2():
    amine2 = stk.StructUnit2.smiles_init(smiles='NCCCN',
                                         functional_groups=['amine'],
                                         name='amine2')
    # Make a second conformer with a distinct geometry.
    amine2.mol.AddConformer(amine2.mol.GetConformer(), True)
    amine2.set_position_from_matrix(
        pos_mat=amine2.mol.GetConformer().GetPositions().T*4,
        conformer=1
    )
    return amine2


def aldehyde2():

    return stk.StructUnit2.smiles_init(smiles='O=CCC=O',
                                       functional_groups=['aldehyde'],
                                       name='aldehyde2')


def polymer(amine2, aldehyde2):
    return stk.Polymer([amine2, aldehyde2],
                       stk.Linear('AB', [0, 0], 3),
                       'polymer')

# gfnxtb = pytest.mark.skipif(
#     all('gfnxtb' not in x for x in sys.argv),
#     reason="Only run when explicitly asked.")
#
#
# @gfnxtb


def test_gfnxtb_properties(tmp_polymer, gfnxtb_path):
    init_dir = os.getcwd()
    # GFNXTB  requires an embedding before working.
    etkdg = stk.ETKDG()
    etkdg.optimize(tmp_polymer)

    # hessian requires optimized structure
    gfnxtb = stk.GFNXTB(gfnxtb_path, output_dir=join(odir, 'gfnxtb_opt'),
                        mem_ulimit=True, opt_level='verytight',
                        num_cores=2)
    gfnxtb.optimize(tmp_polymer)
    energy_calculator = stk.GFNXTBEnergy(gfnxtb_path=gfnxtb_path,
                                         output_dir=join(odir, 'gfnxtb_ey'),
                                         mem_ulimit=True, free=True)
    prop = energy_calculator.energy(tmp_polymer)

    energy = prop['totalenergy']
    freeenergy = prop['totalfreeenergy']
    assert energy != freeenergy
    # check directory moving worked
    assert os.getcwd() == init_dir
    assert len(prop) == 11
    assert isinstance(prop['totalenergy'], float)
    assert isinstance(prop['totalfreeenergy'], float)
    assert isinstance(prop['frequencies'], list)
    assert len(prop['frequencies']) > 0
    for i in prop['frequencies']:
        assert isinstance(i, float)
    assert isinstance(prop['HLGap'], float)
    assert isinstance(prop['FermiLevel'], float)
    assert isinstance(prop['Qdipole'], list)
    assert len(prop['Qdipole']) == 3
    for i in prop['Qdipole']:
        assert isinstance(i, float)
    assert isinstance(prop['fulldipole'], list)
    assert len(prop['fulldipole']) == 4
    for i in prop['fulldipole']:
        assert isinstance(i, float)
    assert isinstance(prop['Qquadrupole'], list)
    assert len(prop['Qquadrupole']) == 6
    for i in prop['Qquadrupole']:
        assert isinstance(i, float)
    assert isinstance(prop['QDIPquadrupole'], list)
    assert len(prop['QDIPquadrupole']) == 6
    for i in prop['QDIPquadrupole']:
        assert isinstance(i, float)
    assert isinstance(prop['fullquadrupole'], list)
    assert len(prop['fullquadrupole']) == 6
    for i in prop['fullquadrupole']:
        assert isinstance(i, float)
    assert prop['occupancies'] is None


def test_valid_solvent():
    gfn_1 = '1'
    gfn_2 = '2'
    valid = 'acetone'  # valid in both versions
    valid1 = 'benzene'  # valid in GFN 1 only
    valid2 = 'dmf'  # valid in GFN 2 only
    invalid = 'andrewtarziawrotethis'
    assert stk.valid_GFNXTB_solvent(solvent=valid, gfn_version=gfn_1)
    assert stk.valid_GFNXTB_solvent(solvent=valid, gfn_version=gfn_2)
    assert stk.valid_GFNXTB_solvent(solvent=valid1, gfn_version=gfn_1)
    assert stk.valid_GFNXTB_solvent(solvent=valid2, gfn_version=gfn_2)
    try:
        stk.valid_GFNXTB_solvent(solvent=valid1, gfn_version=gfn_2)
        assert False
    except stk.GFNXTBInvalidSolventError:
        assert True
    try:
        stk.valid_GFNXTB_solvent(solvent=valid2, gfn_version=gfn_1)
        assert False
    except stk.GFNXTBInvalidSolventError:
        assert True
    try:
        stk.valid_GFNXTB_solvent(solvent=invalid, gfn_version=gfn_1)
        assert False
    except stk.GFNXTBInvalidSolventError:
        assert True
    try:
        stk.valid_GFNXTB_solvent(solvent=invalid, gfn_version=gfn_2)
        assert False
    except stk.GFNXTBInvalidSolventError:
        assert True


def test_gfnxtb_opt(tmp_polymer, gfnxtb_path):
    init_dir = os.getcwd()
    # GFNXTB  requires an embedding before working.
    etkdg = stk.ETKDG()
    etkdg.optimize(tmp_polymer)

    # If the optimization was successful the energy should be lowered.
    energy_calculator = stk.GFNXTBEnergy(gfnxtb_path=gfnxtb_path,
                                         output_dir=join(odir, 'gfnxtb_ey'),
                                         mem_ulimit=True)
    init_prop = energy_calculator.energy(tmp_polymer)
    init_energy = init_prop['totalenergy']

    tmp_polymer.write(join(odir, 'gfnxtb_opt_before.mol'))
    gfnxtb = stk.GFNXTB(gfnxtb_path, output_dir=join(odir, 'gfnxtb_opt'),
                        mem_ulimit=True)
    gfnxtb.optimize(tmp_polymer)
    tmp_polymer.write(join(odir, 'gfnxtb_opt_after.mol'))
    # check directory moving worked
    assert os.getcwd() == init_dir
    properties = energy_calculator.energy(tmp_polymer)
    # optimized structure has lower energy than initial structure
    assert properties['totalenergy'] < init_energy
    energy_calculator = stk.GFNXTBEnergy(gfnxtb_path=gfnxtb_path,
                                         output_dir=join(odir, 'gfnxtb_ey_chrg'),
                                         mem_ulimit=True, charge='-1')
    # energy of structure with formal charge differs from uncharged
    properties_chrg = energy_calculator.energy(tmp_polymer)
    assert properties_chrg['totalenergy'] != init_energy
    energy_calculator = stk.GFNXTBEnergy(gfnxtb_path=gfnxtb_path,
                                         output_dir=join(odir, 'gfnxtb_ey_solv'),
                                         mem_ulimit=True, solvent='h2o')
    # energy of structure in solvent differs from gas phase
    properties_solv = energy_calculator.energy(tmp_polymer)
    assert properties_solv['totalenergy'] != init_energy


def main():
    gfnxtb_path = '/home/atarzia/software/xtb_190418/bin/xtb'
    amine2_ = amine2()
    aldehyde2_ = aldehyde2()
    poly = polymer(amine2_, aldehyde2_)
    test_valid_solvent()
    try:
        test_gfnxtb_properties(tmp_polymer=poly, gfnxtb_path=gfnxtb_path)
    except NotImplementedError:
        print('properties not implemented')
        pass
    try:
        test_gfnxtb_opt(tmp_polymer=poly, gfnxtb_path=gfnxtb_path)
    except NotImplementedError:
        print('opt not implemented')
        pass


if __name__ == '__main__':
    main()
