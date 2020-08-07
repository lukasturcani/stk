import pytest
import numpy as np
import stk

from .case_data import CaseData

bb1 = stk.BuildingBlock('BrCCBr', [stk.BromoFactory()])
bb2 = stk.BuildingBlock('BrCC(CBr)CBr', [stk.BromoFactory()])
topology_graph = stk.cof.PeriodicHoneycomb(
    building_blocks=(bb1, bb2),
    lattice_size=(3, 1, 2),
)


@pytest.fixture(
    params=(
        CaseData(
            periodic_info=topology_graph.get_periodic_info(),
            cell=(
                np.array([109.29499828, 0., 0.]),
                np.array([18.21583305, 31.54982284, 0.]),
                np.array([0., 0., 210.33234855])
            ),
            lengths_and_angles=(
                109.2949982759018,
                36.43086458649658,
                210.3323485478163,
                90.0,
                90.0,
                59.99927221917263,
            )
        ),
    ),
)
def periodic_case(request):

    return request.param
