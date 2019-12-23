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
        assert is_equivalent_atom(a1, a2)

    bonders = it.zip_longest(
        functional_group.get_bonders(),
        clone.get_bonders(),
    )
    for b1, b2 in bonders:
        assert is_equivalent_atom(b1, b2)

    deleters = it.zip_longest(
        functional_group.get_deleters(),
        clone.get_deleters(),
    )
    for d1, d2 in deleters:
        assert is_equivalent_atom(d1, d2)


def is_equivalent_atom(atom1, atom2):
    return (
        atom1 is not atom2
        and atom1.id == atom2.id
        and atom1.__class__ is atom2.__class__
    )


def test_get_atoms(make_functional_group, functional_group_atoms):
    functional_group = make_functional_group(
        atoms=functional_group_atoms.atoms,
        bonders=functional_group_atoms.bonders,
        deleters=functional_group_atoms.deleters,
    )
    atoms = it.zip_longest(
        functional_group.get_atoms(),
        functional_group_atoms.atoms,
    )
    for atom1, atom2 in atoms:
        assert is_equivalent_atom(atom1, atom2)


def test_get_atom_ids(make_functional_group, functional_group_atoms):
    functional_group = make_functional_group(
        atoms=functional_group_atoms.atoms,
        bonders=functional_group_atoms.bonders,
        deleters=functional_group_atoms.deleters,
    )
    atoms = it.zip_longest(
        functional_group.get_atom_ids(),
        functional_group_atoms.atoms,
    )
    for id_, atom in atoms:
        assert id_ == atom.id


def test_get_bonders(make_functional_group, functional_group_atoms):
    functional_group = make_functional_group(
        atoms=functional_group_atoms.atoms,
        bonders=functional_group_atoms.bonders,
        deleters=functional_group_atoms.deleters,
    )
    atoms = it.zip_longest(
        functional_group.get_bonders(),
        functional_group_atoms.bonders,
    )
    for atom1, atom2 in atoms:
        assert is_equivalent_atom(atom1, atom2)


def test_get_bonder_ids(make_functional_group, functional_group_atoms):
    functional_group = make_functional_group(
        atoms=functional_group_atoms.atoms,
        bonders=functional_group_atoms.bonders,
        deleters=functional_group_atoms.deleters,
    )
    atoms = it.zip_longest(
        functional_group.get_bonder_ids(),
        functional_group_atoms.bonders,
    )
    for id_, atom in atoms:
        assert id_ == atom.id


def test_get_deleters(make_functional_group, functional_group_atoms):
    functional_group = make_functional_group(
        atoms=functional_group_atoms.atoms,
        bonders=functional_group_atoms.bonders,
        deleters=functional_group_atoms.deleters,
    )
    atoms = it.zip_longest(
        functional_group.get_deleters(),
        functional_group_atoms.deleters,
    )
    for atom1, atom2 in atoms:
        assert is_equivalent_atom(atom1, atom2)


def test_get_deleter_ids(
    make_functional_group,
    functional_group_atoms,
):
    functional_group = make_functional_group(
        atoms=functional_group_atoms.atoms,
        bonders=functional_group_atoms.bonders,
        deleters=functional_group_atoms.deleters,
    )
    atoms = it.zip_longest(
        functional_group.get_deleter_ids(),
        functional_group_atoms.deleters,
    )
    for id_, atom in atoms:
        assert id_ == atom.id
