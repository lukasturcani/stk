from .utilities import are_equivalent_bonds


def test_get_deleted_bonds(case_data):
    """
    Test that the correct bonds are deleted by a :class:`.Reaction`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the reaction to test and the bonds which
        should be deleted.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_deleted_bonds(
        reaction_result=case_data.reaction.get_result(),
        deleted_bonds=case_data.deleted_bonds,
    )


def _test_get_deleted_bonds(reaction_result, deleted_bonds):
    """
    Test that the correct bonds are deleted by `reaction_result`.

    Parameters
    ----------
    reaction_result : :class:`.ReactionResult`
        The result of a reaction.

    new_bonds : :class:`tuple` of :class:`.Bond`
        The bonds, which should be deleted by `reaction_result`.

    Returns
    -------
    None : :class:`NoneType`

    """

    are_equivalent_bonds(
        reaction_result.get_deleted_bonds(),
        deleted_bonds,
    )
