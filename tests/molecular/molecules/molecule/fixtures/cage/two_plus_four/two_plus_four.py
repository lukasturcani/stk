import pytest
import stk

from ....case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.TwoPlusFour(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles="BrC1=C(Br)[C+]=N1",
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles=(
                                "Br[C+]1C2(Br)[C+]=N[C+]2[C+](Br)[C+]("
                                "Br)[C+2]1"
                            ),
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                ),
            ),
            smiles=(
                "[C+]1=NC2=C1C13[C+]=N[C+]1[C+]1C4=C(N=[C+]4)[C+]4[C+]"
                "2[C+2][C+]2C5=C([C+]=N5)[C+]3[C+2][C+]1C1=C(N=[C+]1)C"
                "21[C+]=N[C+]41"
            ),
            name=name,
        ),
    ),
)
def cage_two_plus_four(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
