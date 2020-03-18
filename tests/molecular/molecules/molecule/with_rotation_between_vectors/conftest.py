import pytest
import stk
import numpy as np


@pytest.fixture(
    params=(
        stk.BuildingBlock('NCCN'),
        stk.BuildingBlock('Brc1ccc(Br)cc1Br', [stk.BromoFactory()]),
        stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(
                    stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
                    stk.BuildingBlock('BrCNCCBr', [stk.BromoFactory()]),
                ),
                repeating_unit='AB',
                num_repeating_units=2,
            ),
        ),
    )
)
def molecule(request):
    """
    A :class:`.Molecule` instance which gets rotated.

    The molecule must have at least 2 atoms for the test to work.

    """

    return request.param.clone()


@pytest.fixture(
    params=[
        np.array([1., 0., 0.]),
        np.array([0., 1., 0.]),
        np.array([0., 0., 1.]),
        np.array([1., 1., 1.]),
    ],
)
def target(request):
    """
    The target vector onto which a molecule is rotated.

    """

    return request.param
