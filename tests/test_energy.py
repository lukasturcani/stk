import stk


def test_rdkit(amine2):
    e1 = amine2.energy.rdkit('uff')
    assert abs(e1 - 37.910347343880076) < 1e-4
    e2 = amine2.energy.rdkit('mmff')
    assert abs(e2 - 7.031976893872189) < 1e-4


def test_formation(polymer, amine2):
    reactant_energy = sum(
        mol.energy.rdkit('mmff')*n
        for mol, n in polymer.bb_counter.items()
    )
    product_energy = (polymer.energy.rdkit('mmff') +
                      amine2.energy.rdkit('mmff'))
    form_energy = reactant_energy - product_energy

    func = stk.FunctionData(name='rdkit',
                            forcefield='mmff')
    calc_form_energy = polymer.energy.formation(
                           func=func,
                           products=[(1, amine2)])

    assert abs(form_energy - calc_form_energy) < 1e-4


def test_pseudoformation(polymer):
    reactant_energy = sum(
            mol.energy.rdkit('mmff')*n
            for mol, n in polymer.bb_counter.items()
    )
    product_energy = polymer.energy.rdkit('mmff')

    pseudoform_energy = reactant_energy - product_energy

    calc_pseudoform_energy = polymer.energy.pseudoformation(
                                stk.FunctionData(name='rdkit',
                                                 forcefield='mmff'))

    assert abs(pseudoform_energy - calc_pseudoform_energy) < 1e-4


def test_logging(amine2):
    amine2.energy.rdkit('uff')
    fd = stk.FunctionData(name='rdkit', forcefield='uff', conformer=-1)
    assert fd in amine2.energy.values
