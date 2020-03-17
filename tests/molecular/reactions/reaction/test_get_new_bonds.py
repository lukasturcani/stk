from .utilities import are_equivalent_bonds


def test_get_new_bonds(case_data):
    """
    Test that correct bonds are added by a :class:`.Reaction`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the reaction to test and the bonds which
        should be addeed.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_new_bonds(
        reaction_result=case_data.reaction.get_result(),
        new_bonds=case_data.new_bonds,
    )


def _test_get_new_bonds(reaction_result, new_bonds):
    """
    Test that correct bonds are found in `reaction_result`.

    Parameters
    ----------
    reaction_result : :class:`.ReactionResult`
        The result of a reaction.

    new_bonds : :class:`tuple` of :class:`.Bond`
        The bonds which should be found in `reaction_result`.

    """

    are_equivalent_bonds(reaction_result.get_new_bonds(), new_bonds)
