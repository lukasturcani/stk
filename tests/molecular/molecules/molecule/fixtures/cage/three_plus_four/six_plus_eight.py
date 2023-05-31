import pytest
import stk

from ....case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.SixPlusEight(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles=(
                                "Br[C+]1[C+2][C+](Br)[C+]2[C+](F)[C+2]"
                                "C2(Br)[C+2]1"
                            ),
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
                "F[C+]1[C+2]C23[C+2][C+]4[C+2][C+]([C+]12)[C+]1[C+2][C"
                "+]2[C+]5[C+2][C+]6[C+2]C7([C+2][C+](F)[C+]57)[C+]5[C+"
                "]7[C+2][C+]([C+]8[C+2][C+]9[C+2]C%10([C+2][C+](F)[C+]"
                "8%10)[C+]8[C+2][C+]3[C+]3[C+]%10[C+2][C+]%11[C+]%12[C"
                "+](F)[C+2]C%12([C+2]%10)[C+]%10[C+2][C+]%12[C+]%13[C+"
                "]%14[C+2][C+]%15[C+]%16[C+2][C+]([C+]%17[C+2][C+]([C+"
                "2]C7%18[C+2][C+](F)[C+]%17%18)C7%17[C+]=N[C+]7[C+]([C"
                "+]([C+2][C+]9%17)C79[C+2][C+](F)[C+]7[C+]([C+2][C+](["
                "C+2]9)C87[C+]=N[C+]37)C%103[C+]=N[C+]%133)C3([C+2][C+"
                "](F)[C+]%153)[C+2]%14)C63[C+]=N[C+]3[C+]%16[C+]3[C+2]"
                "[C+]([C+2]C%126[C+2][C+](F)[C+]36)[C+]2[C+]2N=[C+]C1%"
                "112)C41[C+]=N[C+]51"
            ),
            name=name,
        ),
    ),
)
def cage_six_plus_eight(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
