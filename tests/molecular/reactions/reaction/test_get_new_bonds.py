from .utilities import are_equivalent_bonds


def test_get_new_bonds(test_case):
    _test_get_new_bonds(
        reaction_result=test_case.reaction.get_result(),
        new_bonds=test_case.new_bonds,
    )


def _test_get_new_bonds(reaction_result, new_bonds):
    are_equivalent_bonds(reaction_result.new_bonds, new_bonds)
