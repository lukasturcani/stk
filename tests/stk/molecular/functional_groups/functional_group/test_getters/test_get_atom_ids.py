


def test_get_atom_ids(
    get_functional_group,
    atoms,
    get_bonders,
    get_deleters
):
    functional_group = get_functional_group(
        atoms=atoms,
        bonders=get_bonders(atoms),
        deleters=get_deleters(atoms),
    )
    fg_atoms = it.zip_longest(functional_group.get_atom_ids(), atoms)
    for id_, atom in fg_atoms:
        assert id_ == atom.id
