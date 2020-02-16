from .utilities import are_equivalent_atoms


def test_get_deleted_atoms(test_case):
    _test_get_deleted_atoms(
        reaction_result=test_case.reaction.get_result(),
        deleted_atoms=test_case.deleted_atoms,
    )


def _test_get_deleted_atoms(reaction_result, deleted_atoms):
    are_equivalent_atoms(
        atoms1=reaction_result.deleted_atoms,
        atoms2=deleted_atoms,
    )
