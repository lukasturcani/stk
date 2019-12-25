import itertools as it
from .utilities import is_atom_clone


def test_get_atoms(
    get_functional_group,
    atoms,
    get_bonders,
    get_deleters,
):
    functional_group = get_functional_group(
        atoms=atoms,
        bonders=get_bonders(atoms),
        deleters=get_deleters(atoms),
    )
    fg_atoms = it.zip_longest(functional_group.get_atoms(), atoms)
    for atom1, atom2 in fg_atoms:
        is_atom_clone(atom1, atom2)


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
