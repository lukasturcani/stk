import stk
import pytest


@pytest.mark.parametrize(
    argnames=('molecule', 'num_bonds'),
    argvalues=(
        (stk.BuildingBlock('NCCN'), 11),
    ),
)
def test_get_num_bonds(molecule, num_bonds):
    assert molecule.get_num_bonds() == num_bonds
