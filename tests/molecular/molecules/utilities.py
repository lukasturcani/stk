import numpy as np
import itertools as it


def has_same_structure(molecule1, molecule2):
    assert np.all(np.equal(
        molecule1.get_position_matrix(),
        molecule2.get_position_matrix(),
    ))


def get_displacement_vector(molecule, start_atom, end_atom):
    """
    Get the displacement vector between `start_atom` and `end_atom`.

    """

    position1, position2 = (
        molecule.get_atomic_positions((start_atom, end_atom))
    )
    return position2 - position1


def get_num_atom_ids(molecule, get_atom_ids):
    atom_ids = get_atom_ids(molecule)
    if atom_ids is None:
        return molecule.get_num_atoms()
    else:
        return len(tuple(atom_ids))


def is_equivalent_atom(atom1, atom2):
    assert atom1.get_id() == atom2.get_id()
    assert atom1.get_charge() == atom2.get_charge()
    assert atom1.get_atomic_number() == atom2.get_atomic_number()
    assert atom1.__class__ is atom2.__class__


def is_equivalent_bond(bond1, bond2):
    assert bond1.__class__ is bond2.__class__
    assert bond1.get_order() == bond2.get_order()
    assert bond1.get_periodicity() == bond2.get_periodicity()
    is_equivalent_atom(bond1.get_atom1(), bond2.get_atom1())
    is_equivalent_atom(bond1.get_atom2(), bond2.get_atom2())


def is_equivalent_molecule(molecule1, molecule2):
    atoms = it.zip_longest(
        molecule1.get_atoms(),
        molecule2.get_atoms(),
    )
    for atom1, atom2 in atoms:
        is_equivalent_atom(atom1, atom2)

    bonds = it.zip_longest(
        molecule1.get_bonds(),
        molecule2.get_bonds(),
    )
    for bond1, bond2 in bonds:
        is_equivalent_bond(bond1, bond2)


def normalize_ids(molecule, ids):
    if ids is None:
        return range(molecule.get_num_atoms())
    elif isinstance(ids, int):
        return (ids, )
    else:
        return ids
