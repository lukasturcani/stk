import pytest
import stk

from ....case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.FivePlusTen(
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
                "[C+]1=NC2=C1[C+]1[C+2][C+]3C4=C([C+]=N4)C45[C+]=N[C+]"
                "4[C+]4C6=C(N=[C+]6)[C+]3[C+]3N=[C+]C13C1=C(N=[C+]1)[C"
                "+]1[C+]2[C+2][C+]2C3=C(N=[C+]3)[C+]3[C+2][C+]6C7=C(N="
                "[C+]7)[C+]7[C+2][C+]8C9=C(N=[C+]9)[C+]5[C+2][C+]4C4=C"
                "([C+]=N4)C84[C+]=N[C+]4[C+]7C4=C([C+]=N4)C64[C+]=N[C+"
                "]4[C+]3C3=C([C+]=N3)C23[C+]=N[C+]13"
            ),
            name=name,
        ),
    ),
)
def cage_five_plus_ten(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
