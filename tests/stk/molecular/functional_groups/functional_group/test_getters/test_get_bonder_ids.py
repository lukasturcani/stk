

def test_get_bonder_ids(
    get_functional_group,
    atoms,
    get_bonders,
    get_deleters,
):
    bonders = tuple(get_bonders(atoms))
    functional_group = get_functional_group(
        atoms=atoms,
        bonders=bonders,
        deleters=get_deleters(atoms),
    )
    fg_atoms = it.zip_longest(
        functional_group.get_bonder_ids(),
        bonders,
    )
    for id_, atom in fg_atoms:
        assert id_ == atom.id
