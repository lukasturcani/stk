import stk
import pytest


@pytest.mark.parametrize(
    argnames=('molecule', 'num_atoms'),
    argvalues=(
        (stk.BuildingBlock('NCCN'), 12),
    ),
)
def test_get_num_atoms(molecule, num_atoms):
    assert molecule.get_num_atoms() == num_atoms
