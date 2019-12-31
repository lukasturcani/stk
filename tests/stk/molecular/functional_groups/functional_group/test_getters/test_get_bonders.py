

def test_get_bonders(
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
    fg_atoms = it.zip_longest(functional_group.get_bonders(), bonders)
    for atom1, atom2 in fg_atoms:
        is_atom_clone(atom1, atom2)
