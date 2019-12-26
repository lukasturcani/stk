import itertools as it


def test_clone(molecule):
    clone = molecule.clone()
    atoms = it.zip_longest(clone.get_atoms(), molecule.get_atoms())
    for a1, a2 in atoms:
        is_equivalent_atom(a1, a2)

    bonds = it.zip_longest(clone.get_bonds(), molecule.get_bonds())
    for b1, b2 in bonds:
        is_equivalent_bond(b1, b2)


def is_equivalent_atom(atom1, atom2):
    assert atom1 is not atom2
    assert atom1.id == atom2.id
    assert atom1.charge == atom2.charge
    assert atom1.__class__ is atom2.__class__


def is_equivalent_bond(bond1, bond2):
    assert bond1 is not bond2
    assert bond1.__class__ is bond2.__class__
    assert bond1.order == bond2.order
    assert bond1.periodicity == bond2.periodicity
    is_equivalent_atom(bond1.atom1, bond2.atom1)
    is_equivalent_atom(bond1.atom2, bond2.atom2)
