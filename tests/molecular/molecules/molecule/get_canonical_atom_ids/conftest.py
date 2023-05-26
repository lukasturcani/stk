import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    params=(
        lambda: CaseData(
            molecule=stk.BuildingBlock("[C+2][N+]Br"),
            canonical_atom_ids={0: 0, 1: 2, 2: 1},
        ),
        lambda: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles="Br[C+2][C+2]Br",
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    repeating_unit="A",
                    num_repeating_units=2,
                ),
            ),
            canonical_atom_ids={
                0: 0,
                1: 2,
                2: 4,
                3: 5,
                4: 3,
                5: 1,
            },
        ),
    ),
)
def case_data(request) -> CaseData:
    return request.param()
