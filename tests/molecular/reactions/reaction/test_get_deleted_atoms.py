from .utilities import are_equivalent_atoms


def test_get_deleted_atoms(case_data):
    _test_get_deleted_atoms(
        reaction_result=case_data.reaction.get_result(),
        deleted_atoms=case_data.deleted_atoms,
    )


def _test_get_deleted_atoms(reaction_result, deleted_atoms):
    are_equivalent_atoms(
        atoms1=reaction_result.get_deleted_atoms(),
        atoms2=deleted_atoms,
    )
