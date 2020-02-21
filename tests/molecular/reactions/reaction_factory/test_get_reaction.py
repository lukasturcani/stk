import itertools as it


def test_get_reaction(test_case):
    _test_get_reaction(
        factory=test_case.factory,
        construction_state=test_case.construction_state,
        edges=test_case.edges,
        edge_group=test_case.edge_group,
        functional_groups=test_case.functional_groups,
        reaction_result=test_case.reaction_result,
    )


def _test_get_reaction(
    factory,
    construction_state,
    edges,
    edge_group,
    functional_groups,
    reaction_result,
):
    reaction = factory.get_reaction(
        construction_state=construction_state,
        edges=edges,
        edge_group=edge_group,
        functional_groups=functional_groups,
    )
    is_same_result(reaction.get_result(), reaction_result)


def is_same_result(result1, result2):
    are_same_atoms(result1.new_atoms, result2.new_atoms)
    are_same_bonds(result1.new_bonds, result2.new_bonds)
    are_same_atoms(result1.deleted_atoms, result2.deleted_atoms)


def are_same_atoms(atoms1, atoms2):
    for atom1, atom2 in it.zip_longest(atoms1, atoms2):
        is_same_atom(atom1, atom2)


def is_same_atom(atom1, atom2):
    assert atom1.__class__ is atom2.__class__
    assert atom1.get_id() == atom2.get_id()
    assert atom1.get_charge() == atom2.get_charge()


def are_same_bonds(bonds1, bonds2):
    bonds1 = sorted(bonds1, key=atom_ids)
    bonds2 = sorted(bonds2, key=atom_ids)
    for bond1, bond2 in it.zip_longest(bonds1, bonds2):
        is_same_bond(bond1, bond2)


def atom_ids(bond):
    return bond.get_atom1().get_id(), bond.get_atom2().get_id()


def is_same_bond(bond1, bond2):
    are_same_atoms(
        atoms1=sorted(
            (bond1.get_atom1(), bond1.get_atom2()),
            key=lambda atom: atom.get_id(),
        ),
        atoms2=sorted(
            (bond2.get_atom1(), bond2.get_atom2()),
            key=lambda atom: atom.get_id(),
        )
    )
    assert bond1.get_periodicity() == bond2.get_periodicity()
