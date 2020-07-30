import pytest
import numpy as np
import stk

from .case_data import CaseData

bb1 = stk.BuildingBlock('BrCCBr', [stk.BromoFactory()])
bb2 = stk.BuildingBlock('BrCC(CBr)CBr', [stk.BromoFactory()])


@pytest.fixture(
    params=(
        CaseData(
            topology_graph=stk.cof.Honeycomb(
                building_blocks=(bb1, bb2),
                lattice_size=(3, 1, 2),
                periodic=True,
            ),
            cell=(
                np.array([109.29499828, 0., 0.]),
                np.array([18.21583305, 31.54982284, 0.]),
                np.array([0., 0., 210.33234855])
            ),
        ),
        CaseData(
            topology_graph=stk.cof.Honeycomb(
                building_blocks=(bb1, bb2),
                lattice_size=(3, 3, 1),
                periodic=False,
            ),
            cell=None,
        ),
    ),
)
def periodic_case(request):
    """
    The periodicity of a bond made by a reaction.

    """

    return request.param
