import itertools as it


def is_clone_sequence(atoms, clones, atom_map):
    if atom_map is None:
        atom_map = {}

    for atom, clone in it.zip_longest(atoms, clones):
        is_atom_clone(atom_map.get(atom.id, atom), clone)


def is_atom_clone(atom, clone):
    assert atom is not clone
    assert atom.id == clone.id
    assert atom.charge == clone.charge
    assert atom.__class__ is clone.__class__


def test_clone(functional_group, get_atom_map):
    functional_group.attr = 1
    functional_group._attr = 2
    atom_map = get_atom_map(functional_group)
    clone = functional_group.clone(atom_map)
    assert clone.attr == functional_group.attr
    assert not hasattr(clone, '_attr')
    is_clone_sequence(
        functional_group.get_atoms(),
        clone.get_atoms(),
        atom_map,
    )
    is_clone_sequence(
        functional_group.get_bonders(),
        clone.get_bonders(),
        atom_map,
    )
    is_clone_sequence(
        functional_group.get_deleters(),
        clone.get_deleters(),
        atom_map,
    )


def test_get_atoms(get_functional_group, atoms, bonders, deleters):
    functional_group = get_functional_group(atoms, bonders, deleters)
    fg_atoms = it.zip_longest(functional_group.get_atoms(), atoms)
    for atom1, atom2 in fg_atoms:
        is_atom_clone(atom1, atom2)


def test_get_atom_ids(get_functional_group, atoms, bonders, deleters):
    functional_group = get_functional_group(atoms, bonders, deleters)
    fg_atoms = it.zip_longest(functional_group.get_atom_ids(), atoms)
    for id_, atom in fg_atoms:
        assert id_ == atom.id


def test_get_bonders(get_functional_group, atoms, bonders, deleters):
    functional_group = get_functional_group(atoms, bonders, deleters)
    fg_atoms = it.zip_longest(functional_group.get_bonders(), bonders)
    for atom1, atom2 in fg_atoms:
        is_atom_clone(atom1, atom2)


def test_get_bonder_ids(
    get_functional_group,
    atoms,
    bonders,
    deleters,
):
    functional_group = get_functional_group(atoms, bonders, deleters)
    fg_atoms = it.zip_longest(
        functional_group.get_bonder_ids(),
        bonders,
    )
    for id_, atom in fg_atoms:
        assert id_ == atom.id


def test_get_deleters(get_functional_group, atoms, bonders, deleters):
    functional_group = get_functional_group(atoms, bonders, deleters)
    fg_atoms = it.zip_longest(
        functional_group.get_deleters(),
        deleters,
    )
    for atom1, atom2 in fg_atoms:
        is_atom_clone(atom1, atom2)


def test_get_deleter_ids(
    get_functional_group,
    atoms,
    bonders,
    deleters,
):
    functional_group = get_functional_group(atoms, bonders, deleters)
    fg_atoms = it.zip_longest(
        functional_group.get_deleter_ids(),
        deleters,
    )
    for id_, atom in fg_atoms:
        assert id_ == atom.id
