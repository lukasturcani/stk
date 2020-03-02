import itertools as it


def test_get_bonds(test_case):
    _test_get_bonds(test_case.construction_state, test_case.bonds)


def _test_get_bonds(construction_state, bonds):
    bonds_ = it.zip_longest(
        construction_state.get_bonds(),
        bonds,
    )
    for bond1, bond2 in bonds_:
        assert bond1 is bond2
