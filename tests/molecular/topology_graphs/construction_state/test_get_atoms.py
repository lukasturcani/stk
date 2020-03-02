import itertools as it


def test_get_atoms(test_case):
    _test_get_atoms(test_case.construction_state, test_case.atoms)


def _test_get_atoms(construction_state, atoms):
    atoms_ = it.zip_longest(
        construction_state.get_atoms(),
        atoms,
    )
    for atom1, atom2 in atoms_:
        assert atom1 is atom2
