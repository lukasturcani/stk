from .utilities import are_equivalent_atoms


def test_get_deleted_atoms(test_case):
    _test_get_deleted_atoms(
        reaction=test_case.reaction,
        deleted_atoms=test_case.deleted_atoms,
    )


def _test_get_deleted_atoms(reaction, deleted_atoms):
    are_equivalent_atoms(reaction.get_deleted_atoms(), deleted_atoms)
