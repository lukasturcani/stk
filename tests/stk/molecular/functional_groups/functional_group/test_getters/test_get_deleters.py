


def test_get_deleters(
    get_functional_group,
    atoms,
    get_bonders,
    get_deleters,
):
    deleters = tuple(get_deleters(atoms))
    functional_group = get_functional_group(
        atoms=atoms,
        bonders=get_bonders(atoms),
        deleters=deleters,
    )
    fg_atoms = it.zip_longest(
        functional_group.get_deleters(),
        deleters,
    )
    for atom1, atom2 in fg_atoms:
        is_atom_clone(atom1, atom2)
