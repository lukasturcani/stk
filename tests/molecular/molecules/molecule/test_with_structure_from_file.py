import numpy as np
import os
import pytest
import stk


@pytest.fixture(
    params=[
        'molecule.mol',
        'molecule.xyz',
    ],
)
def path(request, tmpdir):
    return os.path.join(tmpdir, request.param)


def test_with_structure_from_file_0(
    molecule,
    get_position_matrix,
    path,
):
    position_matrix = molecule.get_position_matrix()
    _test_with_structure_from_file_0(
        molecule=molecule,
        get_position_matrix=get_position_matrix,
        path=path,
    )
    # Test immutability.
    assert np.all(np.equal(
        position_matrix,
        molecule.get_position_matrix(),
    ))


def _test_with_structure_from_file_0(
    molecule,
    get_position_matrix,
    path,
):
    position_matrix = get_position_matrix(molecule)
    molecule.with_position_matrix(position_matrix).write(path)
    loaded = molecule.with_structure_from_file(path)
    assert np.allclose(
        a=position_matrix,
        b=loaded.get_position_matrix(),
        atol=1e-4,
    )


@pytest.mark.parametrize(
    argnames=('molecule', 'path'),
    argvalues=(
        (stk.BuildingBlock('NCCN'), 'NCCN.mae'),
    ),
)
def test_with_structure_from_file_1(molecule, datadir, path):
    position_matrix = molecule.get_position_matrix()
    _test_with_structure_from_file_1(molecule, str(datadir / path))
    # Test immutability.
    assert np.all(np.equal(
        position_matrix,
        molecule.get_position_matrix(),
    ))


def _test_with_structure_from_file_1(molecule, path):
    new = molecule.with_structure_from_file(path)
    size_diff = abs(
        molecule.get_maximum_diameter()
        - new.get_maximum_diameter()
    )
    assert size_diff > 1
