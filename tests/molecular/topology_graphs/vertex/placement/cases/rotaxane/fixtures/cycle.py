import pytest
import rdkit.Chem.AllChem as rdkit
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
macrocycle_ids = tuple(max(
    rdkit.GetSymmSSSR(macrocycle.to_rdkit_mol()),
    key=len,
))


def get_plane_normal(building_block):
    macrocycle = max(
        rdkit.GetSymmSSSR(building_block.to_rdkit_mol()),
        key=len,
    )
    normal = building_block.get_plane_normal(macrocycle)
    if np.allclose(normal, [1, 0, 0], atol=1e-13):
        return np.array([1, 0, 0], dtype=np.float64)
    elif np.allclose(normal, [-1, 0, 0], atol=1e-13):
        return np.array([1, 0, 0], dtype=np.float64)

    return normal


@pytest.fixture(
    params=(
        CaseData(
            vertex=vertices._CycleVertex(0, (1, 2, 3), False),
            edges=(),
            building_block=stk.BuildingBlock.init_from_molecule(
                molecule=macrocycle,
            ),
            position=np.array([1, 2, 3], dtype=np.float64),
            alignment_tests={
                get_plane_normal:
                    np.array([1, 0, 0], dtype=np.float64),
            },
            functional_group_edges={},
            position_ids=macrocycle_ids,
        ),
        CaseData(
            vertex=vertices._CycleVertex(0, (1, 2, 3), True),
            edges=(),
            building_block=stk.BuildingBlock.init_from_molecule(
                molecule=macrocycle,
            ),
            position=np.array([1, 2, 3], dtype=np.float64),
            alignment_tests={
                # If correctly aligned, get_plane_normal returns the
                # same vector, regardless of flip status.
                get_plane_normal:
                    np.array([1, 0, 0], dtype=np.float64),
            },
            functional_group_edges={},
            position_ids=macrocycle_ids,
        ),
    ),
)
def cycle(request):
    return request.param
