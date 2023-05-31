import pytest
import stk

from ....case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.TwoPlusTwo(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles=(
                                "Br[C+]1[C+2][C+](Br)[C+]2[C+](F)[C+2]"
                                "C2(Br)[C+2]1"
                            ),
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                ),
            ),
            smiles=(
                "F[C+]1[C+2]C23[C+2][C+]4[C+2][C+]([C+]5[C+2][C+]6[C+]"
                "7[C+](F)[C+2]C7([C+2]5)[C+]5[C+2][C+]2[C+2]C2([C+2][C"
                "+](F)[C+]52)[C+]2[C+2][C+]4[C+2]C64[C+2][C+](F)[C+]24"
                ")[C+]13"
            ),
            name=name,
        ),
    ),
)
def cage_two_plus_two(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
