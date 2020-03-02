import itertools as it


def test_get_atom_infos(test_case):
    _test_get_atom_infos(
        construction_state=test_case.construction_state,
        atom_infos=test_case.atom_infos,
    )


def _test_get_atom_infos(construction_state, atom_infos):
    infos = it.zip_longest(
        construction_state.get_atom_infos(),
        atom_infos,
    )
    for atom_info1, atom_info2 in infos:
        atom_info1 is atom_info2
