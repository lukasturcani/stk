import pytest
import numpy as np
import stk

from .case_data import CaseData

bb1 = stk.BuildingBlock('BrCCBr', [stk.BromoFactory()])
bb2 = stk.BuildingBlock('BrCC(CBr)CBr', [stk.BromoFactory()])


@pytest.fixture(
    params=(
        CaseData(
            topology_graph=stk.cof.PeriodicHoneycomb(
                building_blocks=(bb1, bb2),
                lattice_size=(3, 1, 2),
            ),
            cell=(
                np.array([109.29499828, 0., 0.]),
                np.array([18.21583305, 31.54982284, 0.]),
                np.array([0., 0., 210.33234855])
            ),
        ),
        CaseData(
            topology_graph=stk.cof.PeriodicHoneycomb(
                building_blocks=(bb1, bb2),
                lattice_size=(3, 1, 2),
                optimizer=stk.PeriodicCollapser(),
            ),
            cell=(
                np.array([109.29499828, 0., 0.]),
                np.array([18.21583305, 31.54982284, 0.]),
                np.array([0., 0., 210.33234855])
            ),
        ),
    ),
)
def unscaled_periodic_case(request):

    return request.param


@pytest.fixture(
    params=(
        CaseData(
            topology_graph=stk.cof.PeriodicHoneycomb(
                building_blocks=(bb1, bb2),
                lattice_size=(3, 1, 2),
            ),
            cell=(
                np.array([109.29499828, 0., 0.]),
                np.array([18.21583305, 31.54982284, 0.]),
                np.array([0., 0., 210.33234855])
            ),
        ),
        CaseData(
            topology_graph=stk.cof.PeriodicHoneycomb(
                building_blocks=(bb1, bb2),
                lattice_size=(3, 1, 2),
                optimizer=stk.PeriodicCollapser(),
            ),
            cell=(
                np.array([40.32953729, 0., 0.]),
                np.array([7.88334589, 13.65395508, 0.]),
                np.array([0., 0., 76.68756542])
            ),
        ),
    ),
)
def scaled_periodic_case(request):

    return request.param
