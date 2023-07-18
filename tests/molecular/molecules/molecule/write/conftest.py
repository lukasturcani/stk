import os
import pathlib

import pytest


@pytest.fixture(
    params=[
        "molecule.mol",
        "molecule.xyz",
        pathlib.Path("molecule-path.mol"),
        pathlib.Path("molecule-path.xyz"),
    ],
)
def path(request, tmpdir):
    return os.path.join(tmpdir, request.param)
