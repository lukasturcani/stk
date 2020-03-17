import itertools as it


def are_equivalent_atoms(atoms1, atoms2):
    for atom1, atom2 in it.zip_longest(atoms1, atoms2):
        is_equivalent_atom(atom1, atom2)


def is_equivalent_atom(atom1, atom2):
    assert atom1.get_id() == atom2.get_id()
    assert atom1.__class__ is atom2.__class__
    assert atom1.get_atomic_number() == atom2.get_atomic_number()
    assert atom1.get_charge() == atom2.get_charge()


def are_equivalent_bonds(bonds1, bonds2):
    for bond1, bond2 in it.zip_longest(bonds1, bonds2):
        is_equivalent_bond(bond1, bond2)


def is_equivalent_bond(bond1, bond2):
    assert bond1.__class__ is bond2.__class__
    assert bond1.get_order() == bond2.get_order()
    assert bond1.get_periodicity() == bond2.get_periodicity()
    is_equivalent_atom(bond1.get_atom1(), bond2.get_atom1())
    is_equivalent_atom(bond1.get_atom2(), bond2.get_atom2())
