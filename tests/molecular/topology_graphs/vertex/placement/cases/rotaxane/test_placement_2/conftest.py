import pytest
import rdkit.Chem.AllChem as rdkit
import numpy as np
import stk

from .case_data import CaseData

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
            vertex1=vertices._CycleVertex(0, (1, 2, 3), False),
            vertex2=vertices._CycleVertex(0, (1, 2, 3), True),
            building_block=stk.BuildingBlock.init_from_molecule(
                molecule=macrocycle,
            ),
            atom_ids=macrocycle_ids,
        ),
    ),
)
def case_data(request):
    return request.param
