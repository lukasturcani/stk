import pytest
import numpy as np


@pytest.fixture(
    params=[
        [0, 0, 0],
        [10, 20, 30],
        [-10, 20, -30],
        [0.5, 10, -0.921],
    ],
)
def displacement(request):
    return list(request.param)


def test_with_displacement(molecule, displacement):
    new = molecule.with_displacement(displacement)
    assert np.allclose(
        a=molecule.get_position_matrix()+displacement,
        b=new.get_position_matrix(),
        atol=1e-32,
    )
