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
            x_vector=np.array([109.29499828, 0., 0.]),
            y_vector=np.array([18.21583305, 31.54982284, 0.]),
            z_vector=np.array([0., 0., 210.33234855]),
            a=109.2949982759018,
            b=36.43086458649658,
            c=210.3323485478163,
            alpha=90.0,
            beta=90.0,
            gamma=59.99927221917263,
        ),
    ),
)
def periodic_case(request):

    return request.param
