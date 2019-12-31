

def test_get_deleter_ids(
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
        functional_group.get_deleter_ids(),
        deleters,
    )
    for id_, atom in fg_atoms:
        assert id_ == atom.id
