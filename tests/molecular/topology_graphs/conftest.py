import numpy as np
import pytest

import stk

from .case_data import CaseData

bb1 = stk.BuildingBlock('BrCCBr', [stk.BromoFactory()])
bb2 = stk.BuildingBlock('BrCC(CBr)CBr', [stk.BromoFactory()])


@pytest.fixture(
    params=(
        lambda: CaseData(
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
        lambda: CaseData(
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
def unscaled_periodic_case(request) -> CaseData:
    """
    A :class:`.CaseData` instance.

    The instance produced by this fixture is unscaled because the
    expected cell constants are not scaled, even if the molecule
    is used with an optimizer. This is used for testing
    :meth:`.TopologyGraph.get_periodic_info`, which does not scale
    the cell constants.

    """

    return request.param()


@pytest.fixture(
    params=(
        lambda: CaseData(
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
        lambda: CaseData(
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
def scaled_periodic_case(request) -> CaseData:
    """
    A :class:`.CaseData` instance.

    The instance produced by this fixture is scaled because the
    expected cell constants are scaled, assuming the collapser
    performs scaling on the unit cell. This is used for testing
    :meth:`.PeriodicConstructionResult.get_periodic_info`, which
    returns scaled cell constants.

    """

    return request.param()
