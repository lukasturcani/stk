import pytest
import numpy as np
import stk

from ....case_data import CaseData

vertices = stk.molecular.topology_graphs.rotaxane.vertices


macrocycle = stk.ConstructedMolecule(
    topology_graph=stk.macrocycle.Macrocycle(
        building_blocks=(
            stk.BuildingBlock('BrCCCBr', [stk.BromoFactory()]),
        ),
        repeating_unit='A',
        num_repeating_units=5,
    ),
)


@pytest.fixture(
    params=(
        CaseData(
            vertex=vertices._CycleVertex(0, (1, 2, 3), True),
            edges=(),
            building_block=stk.BuildingBlock.init_from_molecule(
                molecule=macrocycle,
            ),
            position=np.array([1, 2, 3], dtype=np.float64),
            alignment_tests={},
            functional_group_edges={},
        ),
        CaseData(
            vertex=vertices._CycleVertex(0, (1, 2, 3), False),
            edges=(),
            building_block=stk.BuildingBlock.init_from_molecule(
                molecule=macrocycle,
            ),
            position=np.array([1, 2, 3], dtype=np.float64),
            alignment_tests={},
            functional_group_edges={},
        ),
    ),
)
def axle(request):
    return request.param
