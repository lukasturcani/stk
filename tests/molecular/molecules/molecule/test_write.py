import numpy as np
import os
import pytest


@pytest.fixture(
    params=[
        'molecule.mol',
        'molecule.xyz',
    ],
)
def path(request, tmpdir):
    return os.path.join(tmpdir, request.param)


def test_write(molecule, get_position_matrix, path,):
    position_matrix = get_position_matrix(molecule)
    molecule.with_position_matrix(position_matrix).write(path)
    loaded = molecule.with_structure_from_file(path)
    assert np.allclose(
        a=position_matrix,
        b=loaded.get_position_matrix(),
        atol=1e-4,
    )
