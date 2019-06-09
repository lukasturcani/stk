import stk
from os.path import join
import os

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


def test_gfnxtb_energy(tmp_polymer, gfnxtb_path):
    raise NotImplementedError


def test_gfnxtb_properties(tmp_polymer, gfnxtb_path):
    raise NotImplementedError


def test_gfnxtb_ehess(tmp_polymer, gfnxtb_path):
    raise NotImplementedError


def test_gfnxtb_optsolv(tmp_polymer, gfnxtb_path):
    raise NotImplementedError


def test_gfnxtb_optchrg(tmp_polymer, gfnxtb_path):
    raise NotImplementedError


def test_gfnxtb_opt(tmp_polymer, gfnxtb_path):
    # GFNXTB  requires an embedding before working.
    etkdg = stk.ETKDG()
    etkdg.optimize(tmp_polymer)

    # If the optimization was successful the energy should be lowered.
    energy_calculator = stk.GFNXTBEnergy(gfnxtb_path=gfnxtb_path)
    init_energy = energy_calculator.energy(tmp_polymer)

    print('doing GFN optimization')

    tmp_polymer.write(join(odir, 'gfnxtb_opt_before.mol'))
    gfnxtb = stk.GFNXTB(gfnxtb_path, output_dir=join(odir, 'gfnxtb_opt'),
                        mem_ulimit=True)
    gfnxtb.optimize(tmp_polymer)
    tmp_polymer.write(join(odir, 'gfnxtb_opt_after.mol'))

    print('done GFN optimization')
    assert energy_calculator.energy(tmp_polymer) < init_energy


def main():
    gfnxtb_path = '/home/atarzia/software/xtb_190418/bin/xtb'
    amine2_ = amine2()
    aldehyde2_ = aldehyde2()
    poly = polymer(amine2_, aldehyde2_)
    try:
        test_gfnxtb_energy(tmp_polymer=poly, gfnxtb_path=gfnxtb_path)
    except NotImplementedError:
        print('energy not implemented')
        pass
    try:
        test_gfnxtb_properties(tmp_polymer=poly, gfnxtb_path=gfnxtb_path)
    except NotImplementedError:
        print('properties not implemented')
        pass
    try:
        test_gfnxtb_ehess(tmp_polymer=poly, gfnxtb_path=gfnxtb_path)
    except NotImplementedError:
        print('ehess not implemented')
        pass
    try:
        test_gfnxtb_opt(tmp_polymer=poly, gfnxtb_path=gfnxtb_path)
    except NotImplementedError:
        print('opt not implemented')
        pass
    try:
        test_gfnxtb_optsolv(tmp_polymer=poly, gfnxtb_path=gfnxtb_path)
    except NotImplementedError:
        print('optsolv not implemented')
        pass
    try:
        test_gfnxtb_optchrg(tmp_polymer=poly, gfnxtb_path=gfnxtb_path)
    except NotImplementedError:
        print('ocharg not implemented')
        pass


if __name__ == '__main__':
    main()
