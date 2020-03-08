import stk
import numpy as np
import pytest

from ..utilities import has_same_structure, get_displacement_vector


@pytest.fixture(
    params=(
        stk.BuildingBlock('NCCN'),
        stk.BuildingBlock('Brc1ccc(Br)cc1Br', [stk.BromoFactory()]),
        stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(
                    stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
                    stk.BuildingBlock(
                        smiles='BrCNCCBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                repeating_unit='AB',
                num_repeating_units=2,
            ),
        ),
    ),
)
def molecule(request):
    """
    A :class:`.Molecule` instance with at least 3 atoms.

    """

    return request.param.clone()


def test_with_rotation_to_minimize_angle(molecule):
    # Use to check that immutability is not violated.
    clone = molecule.clone()
    _test_with_rotation_to_minimize_angle(molecule)
    has_same_structure(molecule, clone)


def _test_with_rotation_to_minimize_angle(molecule):
    start = get_displacement_vector(molecule, 0, 1)
    target = get_displacement_vector(molecule, 0, 2)
    new = molecule.with_rotation_to_minimize_angle(
        start=start,
        target=target,
        axis=stk.normalize_vector(np.cross(start, target)),
        origin=next(molecule.get_atomic_positions((0, ))),
    )

    result = get_displacement_vector(new, 0, 1)
    assert np.allclose(
        a=stk.normalize_vector(result),
        b=stk.normalize_vector(target),
        atol=1e-12,
    )
    assert abs(np.linalg.norm(start) - np.linalg.norm(result)) < 1e-15
