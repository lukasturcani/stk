import itertools as it

import numpy as np


def test_get_reaction(case_data):
    """
    Test :class:`.ReactionFactory.get_reaction`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        The test case. Holds the factory, and the expected reaction
        result.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_reaction(
        factory=case_data.factory,
        construction_state=case_data.construction_state,
        edge_group=case_data.edge_group,
        reaction_result=case_data.reaction_result,
    )


def _test_get_reaction(
    factory,
    construction_state,
    edge_group,
    reaction_result,
):
    """
    Test :class:`.ReactionFactory.get_reaction`.

    Parameters
    ----------
    factory : :class:`.ReactionFactory`
        The factory to test.

    construction_state : :class:`.ConstructionState`
        The construction state to pass to the factory.

    edge_group : :class:`.EdgeGroup`
        The edge group to pass to the factory.

    reaction_result : :class:`.ReactionResult`
        The expected result of the reaction returned by
        `factory`.

    Returns
    -------
    None : :class:`NoneType`

    """

    reaction = factory.get_reaction(
        construction_state=construction_state,
        edge_group=edge_group,
    )
    is_same_result(reaction.get_result(), reaction_result)


def is_same_result(result1, result2):
    are_same_bonds(result1.get_new_bonds(), result2.get_new_bonds())
    are_same_new_atoms(
        new_atoms1=result1.get_new_atoms(),
        new_atoms2=result2.get_new_atoms(),
    )
    are_same_atoms(
        atoms1=result1.get_deleted_atoms(),
        atoms2=result2.get_deleted_atoms(),
    )
    are_same_bonds(
        result1.get_deleted_bonds(),
        result2.get_deleted_bonds(),
    )


def are_same_new_atoms(new_atoms1, new_atoms2):
    for (atom1, position1), (atom2, position2) in it.zip_longest(
        new_atoms1,
        new_atoms2,
    ):
        is_same_atom(atom1, atom2)
        assert np.allclose(position1, position2, atol=1e-13)


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
        ),
    )
    assert bond1.get_periodicity() == bond2.get_periodicity()
