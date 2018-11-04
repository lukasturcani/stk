import stk


def test_rdkit(mol):
    e1 = mol.energy.rdkit('uff')
    assert abs(e1 - 59.003470735439585) < 1e-4
    e2 = mol.energy.rdkit('mmff')
    assert abs(e2 - 65.17285719079524) < 1e-4


def test_formation(polymer, mol):
    reactant_energy = sum(mol.energy.rdkit('uff') for
                          mol in polymer.building_blocks)
    product_energy = (polymer.energy.rdkit('uff') +
                      mol.energy.rdkit('uff'))
    form_energy = reactant_energy - product_energy

    func = stk.FunctionData(name='rdkit',
                            forcefield='uff')
    calc_form_energy = polymer.energy.formation(
                           func=func,
                           products=[(1, mol)])

    assert abs(form_energy - calc_form_energy) < 1e-4


def test_pseudoformation(polymer):
    reactant_energy = sum(mol.energy.rdkit('uff') for
                          mol in polymer.building_blocks)
    product_energy = polymer.energy.rdkit('uff')
    pseudoform_energy = reactant_energy - product_energy

    calc_pseudoform_energy = polymer.energy.pseudoformation(
                                stk.FunctionData(name='rdkit',
                                                 forcefield='uff'))

    assert abs(pseudoform_energy - calc_pseudoform_energy) < 1e-4


def test_logging(mol):
    mol.energy.rdkit('uff')
    fd = stk.FunctionData(name='rdkit', forcefield='uff', conformer=-1)
    assert fd in mol.energy.values
