import itertools as it


def test_clone(functional_group, make_atom_map):
    functional_group.attr = 1
    functional_group._attr = 2
    atom_map = make_atom_map(functional_group)
    clone = functional_group.clone(atom_map)
    assert clone.attr == functional_group.attr
    assert not hasattr(clone, '_attr')

    atoms = it.zip_longest(
        functional_group.get_atoms(),
        clone.get_atoms(),
    )
    for a1, a2 in atoms:
        _test_cloned_atom(a1, a2, atom_map)

    bonders = it.zip_longest(
        functional_group.get_bonders(),
        clone.get_bonders(),
    )
    for b1, b2 in bonders:
        _test_cloned_atom(b1, b2, atom_map)

    deleters = it.zip_longest(
        functional_group.get_deleters(),
        clone.get_deleters(),
    )
    for d1, d2 in deleters:
        _test_cloned_atom(d1, d2, atom_map)


def _test_cloned_atom(atom, clone, atom_map):
    if atom_map is None or atom.id not in atom_map:
        assert atom is not clone
        assert atom.id == clone.id
        assert atom.__class__ is clone.__class__
    else:
        assert atom_map[atom.id] is not clone
        assert atom_map[atom.id].id == clone.id
        assert atom_map[atom.id].__class__ is clone.__class__


def _test_atom(atom1, atom2):
    assert atom1 is not atom2
    assert atom1.id == atom2.id
    assert atom1.__class__ is atom2.__class__


def test_get_atoms(make_functional_group):
    fg_data = make_functional_group()
    atoms = it.zip_longest(
        fg_data.functional_group.get_atoms(),
        fg_data.atoms,
    )
    for atom1, atom2 in atoms:
        _test_atom(atom1, atom2)


def test_get_atom_ids(make_functional_group):
    fg_data = make_functional_group()
    atoms = it.zip_longest(
        fg_data.functional_group.get_atom_ids(),
        fg_data.atoms,
    )
    for id_, atom in atoms:
        assert id_ == atom.id


def test_get_bonders(make_functional_group):
    fg_data = make_functional_group()
    atoms = it.zip_longest(
        fg_data.functional_group.get_bonders(),
        fg_data.bonders,
    )
    for atom1, atom2 in atoms:
        _test_atom(atom1, atom2)


def test_get_bonder_ids(make_functional_group):
    fg_data = make_functional_group()
    atoms = it.zip_longest(
        fg_data.functional_group.get_bonder_ids(),
        fg_data.bonders,
    )
    for id_, atom in atoms:
        assert id_ == atom.id


def test_get_deleters(make_functional_group):
    fg_data = make_functional_group()
    atoms = it.zip_longest(
        fg_data.functional_group.get_deleters(),
        fg_data.deleters,
    )
    for atom1, atom2 in atoms:
        _test_atom(atom1, atom2)


def test_get_deleter_ids(make_functional_group):
    fg_data = make_functional_group()
    atoms = it.zip_longest(
        fg_data.functional_group.get_deleter_ids(),
        fg_data.deleters,
    )
    for id_, atom in atoms:
        assert id_ == atom.id
