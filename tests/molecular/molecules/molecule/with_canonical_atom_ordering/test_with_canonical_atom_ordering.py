from ...utilities import is_clone


def test_with_canonical_atom_ordering(case_data):
    _test_with_canonical_atom_ordering(
        molecule=case_data.molecule,
        result=case_data.result,
    )


def _test_with_canonical_atom_ordering(molecule, result):
    is_clone(molecule.with_canonical_atom_ordering(), result)
