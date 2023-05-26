import pytest
import stk

from ....case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.TwoPlusThree(
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
                "F[C+]1[C+2]C23[C+2][C+]4[C+2][C+](C5=C(N=[C+]5)C56[C+"
                "2][C+]([C+2][C+](C7=C2[C+]=N7)[C+]5[C+](F)[C+2]6)C2=C"
                "4[C+]=N2)[C+]13"
            ),
            name=name,
        ),
    ),
)
def cage_two_plus_three(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
