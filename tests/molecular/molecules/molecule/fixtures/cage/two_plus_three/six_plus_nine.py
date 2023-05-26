import pytest
import stk

from ....case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.SixPlusNine(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles="BrC1=C(Br)[C+]=N1",
                            functional_groups=[stk.BromoFactory()],
                        ),
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
                "F[C+]1[C+2]C23[C+2][C+]4[C+2][C+](C5=C(N=[C+]5)[C+]5["
                "C+2][C+]6[C+2]C7([C+2][C+](F)[C+]57)C5=C(N=[C+]5)[C+]"
                "5[C+2][C+]([C+2]C7([C+2][C+](F)[C+]57)C5=C4[C+]=N5)C4"
                "=C(N=[C+]4)[C+]4[C+2][C+]5C7=C([C+]=N7)[C+]7[C+2][C+]"
                "(C8=C2[C+]=N8)[C+]2[C+](F)[C+2]C2([C+2]7)C2=C(N=[C+]2"
                ")[C+]2[C+2][C+]([C+2]C7([C+2][C+](F)[C+]27)C2=C6[C+]="
                "N2)C2=C(N=[C+]2)C2([C+2]4)[C+2][C+](F)[C+]52)[C+]13"
            ),
            name=name,
        ),
    ),
)
def cage_six_plus_nine(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
