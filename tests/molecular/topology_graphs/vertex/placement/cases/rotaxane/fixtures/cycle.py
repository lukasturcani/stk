import numpy as np
import pytest
import rdkit.Chem.AllChem as rdkit
import stk

from ....case_data import CaseData


def _get_macrocycle() -> stk.BuildingBlock:
    macrocycle = stk.ConstructedMolecule(
        topology_graph=stk.macrocycle.Macrocycle(
            building_blocks=(
                stk.BuildingBlock("BrCCCBr", [stk.BromoFactory()]),
            ),
            repeating_unit="A",
            num_repeating_units=5,
        ),
    )
    return stk.BuildingBlock.init_from_molecule(macrocycle)


def _get_case_data_1() -> CaseData:
    macrocycle = _get_macrocycle()
    macrocycle_ids = tuple(
        max(
            rdkit.GetSymmSSSR(macrocycle.to_rdkit_mol()),
            # mypy throws an incorrect type error here, I believe it has
            # issues with generic functions sometimes.
            key=len,  # type: ignore
        )
    )
    return CaseData(
        vertex=stk.rotaxane.CycleVertex(0, (1, 2, 3), False),
        edges=(),
        building_block=_get_macrocycle(),
        position=np.array([1, 2, 3], dtype=np.float64),
        alignment_tests={
            get_plane_normal: np.array([1, 0, 0], dtype=np.float64),
        },
        functional_group_edges={},
        position_ids=macrocycle_ids,
    )


def _get_case_data_2() -> CaseData:
    macrocycle = _get_macrocycle()
    macrocycle_ids = tuple(
        max(
            rdkit.GetSymmSSSR(macrocycle.to_rdkit_mol()),
            # mypy throws an incorrect type error here, I believe it has
            # issues with generic functions sometimes.
            key=len,  # type: ignore
        )
    )
    return CaseData(
        vertex=stk.rotaxane.CycleVertex(0, (1, 2, 3), True),
        edges=(),
        building_block=_get_macrocycle(),
        position=np.array([1, 2, 3], dtype=np.float64),
        alignment_tests={
            get_plane_normal: np.array([-1, 0, 0], dtype=np.float64),
        },
        functional_group_edges={},
        position_ids=macrocycle_ids,
    )


def get_plane_normal(building_block):
    macrocycle = max(
        rdkit.GetSymmSSSR(building_block.to_rdkit_mol()),
        key=len,
    )
    centroid = building_block.get_centroid()
    (atom1_position,) = building_block.get_atomic_positions(0)
    normal = -stk.get_acute_vector(
        reference=atom1_position - centroid,
        vector=building_block.get_plane_normal(macrocycle),
    )
    if np.allclose(normal, [1, 0, 0], atol=1e-13):
        return np.array([1, 0, 0], dtype=np.float64)
    elif np.allclose(normal, [-1, 0, 0], atol=1e-13):
        return np.array([-1, 0, 0], dtype=np.float64)

    return normal


@pytest.fixture(
    scope="session",
    params=(
        _get_case_data_1,
        _get_case_data_2,
    ),
)
def cycle(request) -> CaseData:
    return request.param()
