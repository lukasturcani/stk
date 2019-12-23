def _test_atom(atom, clone_atom, atom_map):
    if atom_map is None or atom.id not in atom_map:
        assert atom is not clone_atom
        assert atom.id == clone_atom.id
        assert atom.__class__ is clone_atom.__class__
    else:
        assert clone_atom is atom_map[atom.id]


def test_clone(bond, get_atom_map):
    atom_map = get_atom_map(bond)
    bond.attr = 1
    bond._attr = 2
    clone = bond.clone(atom_map)
    _test_atom(bond.atom1, clone.atom1, atom_map)
    _test_atom(bond.atom2, clone.atom2, atom_map)
    assert bond.periodicity == clone.periodicity
    assert bond.attr == clone.attr
    assert not hasattr(clone, '_attr')


def test_is_periodic(bond):
    is_periodic = any(x != 0 for x in bond.periodicity)
    assert is_periodic == bond.is_periodic()
