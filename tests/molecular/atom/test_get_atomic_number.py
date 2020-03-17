def test_get_atomic_number(case_data):
    _test_get_atomic_number(case_data.atom, case_data.atomic_number)


def _test_get_atomic_number(atom, atomic_number):
    assert atom.get_atomic_number() == atomic_number
