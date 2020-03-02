import itertools as it


def test_get_bond_infos(test_case):
    _test_get_bond_infos(
        construction_state=test_case.construction_state,
        bond_infos=test_case.bond_infos,
    )


def _test_get_bond_infos(construction_state, bond_infos):
    infos = it.zip_longest(
        construction_state.get_bond_infos(),
        bond_infos,
    )
    for bond_info1, bond_info2 in infos:
        assert bond_info1 is bond_info2
