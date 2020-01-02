from .utilities import are_equivalent_bonds


def test_get_new_bonds(test_case):
    _test_get_new_bonds(
        reaction=test_case.reaction,
        new_bonds=test_case.new_bonds,
    )


def _test_get_new_bonds(reaction, new_bonds):
    are_equivalent_bonds(reaction.get_new_bonds(), new_bonds)
