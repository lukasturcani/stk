import itertools as it
import stk
import pytest


@pytest.mark.parametrize(
    argnames=('molecule', 'atoms'),
    argvalues=(
        (
            stk.BuildingBlock('NCCN'),
            (
                stk.N(0),
                stk.C(1),
                stk.C(2),
                stk.N(3),
                stk.H(4),
                stk.H(5),
                stk.H(6),
                stk.H(7),
                stk.H(8),
                stk.H(9),
                stk.H(10),
                stk.H(11),
            ),
        ),
    ),
)
def test_get_atoms(molecule, atoms, get_atom_ids):
    atom_ids = get_atom_ids(molecule)
    if atom_ids is None:
        atom_ids = range(molecule.get_num_atoms())

    results = it.zip_longest(
        atom_ids,
        molecule.get_atoms(get_atom_ids(molecule)),
    )
    for atom_id, atom in results:
        is_equivalent_atom(atoms[atom_id], atom)


def is_equivalent_atom(atom1, atom2):
    assert atom1 is not atom2
    assert atom1.id == atom2.id
    assert atom1.charge == atom2.charge
    assert atom1.__class__ is atom2.__class__
