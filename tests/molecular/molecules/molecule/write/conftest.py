import pytest
import os


@pytest.fixture(
    params=[
        'molecule.mol',
        'molecule.xyz',
    ],
)
def path(request, tmpdir):
    return os.path.join(tmpdir, request.param)
