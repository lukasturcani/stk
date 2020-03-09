from .utilities import are_equivalent_bonds


def test_get_new_bonds(case_data):
    _test_get_new_bonds(
        reaction_result=case_data.reaction.get_result(),
        new_bonds=case_data.new_bonds,
    )


def _test_get_new_bonds(reaction_result, new_bonds):
    are_equivalent_bonds(reaction_result.get_new_bonds(), new_bonds)
