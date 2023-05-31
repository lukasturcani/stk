import pytest
import stk

from ....case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.FourPlusFour(
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
                "7[C+](F)[C+2]C7([C+2]5)[C+]5[C+2][C+]7[C+2]C8([C+2][C"
                "+](F)[C+]58)[C+]5[C+2][C+]8[C+2]C9([C+2][C+](F)[C+]59"
                ")[C+]5[C+2][C+]([C+2]C69[C+2][C+](F)[C+]59)C56[C+2][C"
                "+](F)[C+]5[C+]([C+2][C+]4[C+2]6)[C+]4[C+2][C+]([C+]5["
                "C+](F)[C+2]C85[C+2]4)C45[C+2][C+](F)[C+]4[C+]7[C+2][C"
                "+]2[C+2]5)[C+]13"
            ),
            name=name,
        ),
    ),
)
def cage_four_plus_four(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
