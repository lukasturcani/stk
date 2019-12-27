import itertools as it
import stk
import pytest

from ..utilities import is_equivalent_atom, sanitize_ids


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
        sanitize_ids(molecule, atom_ids),
        molecule.get_atoms(get_atom_ids(molecule)),
    )
    for atom_id, atom in results:
        is_equivalent_atom(atoms[atom_id], atom)
