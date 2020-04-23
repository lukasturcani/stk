from .utilities import are_equivalent_atoms


def test_get_deleted_atoms(case_data):
    """
    Test correction deletion of atoms by a :class:`.Reaction`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Includes the reaction and atoms which should have
        been deleted.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_deleted_atoms(
        reaction_result=case_data.reaction.get_result(),
        deleted_atoms=case_data.deleted_atoms,
    )


def _test_get_deleted_atoms(reaction_result, deleted_atoms):
    """
    Test that the correct atoms are deleted by `reaction_result`.

    Parameters
    ----------
    reaction_result : :class:`.ReactionResult`
        The result of a reaction.

    deleted_atoms : :class:`tuple` of :class:`.Atom`
        The atoms which should be deleted by the reaction.

    Returns
    -------
    None : :class:`NoneType`

    """

    are_equivalent_atoms(
        atoms1=reaction_result.get_deleted_atoms(),
        atoms2=deleted_atoms,
    )
