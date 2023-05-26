import pytest
import stk

from ....case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.TwentyPlusThirty(
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
                "2][C+]7[C+2][C+](C8=C([C+]=N8)C89[C+2][C+]%10[C+2][C+"
                "](C%11=C(N=[C+]%11)[C+]%11[C+2][C+]%12[C+2]C%13([C+2]"
                "[C+](F)[C+]%11%13)C%11=C([C+]=N%11)C%11%13[C+2][C+]%1"
                "4[C+2][C+](C%15=C(N=[C+]%15)C%15%16[C+2][C+]%17[C+2]["
                "C+](C%18=C([C+]=N%18)C%18%19[C+2][C+]([C+2][C+](C%20="
                "C(N=[C+]%20)[C+]%20[C+2][C+](C%21=C2[C+]=N%21)[C+]2[C"
                "+](F)[C+2]C2([C+2]%20)C2=C%12[C+]=N2)[C+]%18[C+](F)[C"
                "+2]%19)C2=C(N=[C+]2)C2%12[C+2][C+]%18[C+2][C+](C%19=C"
                "([C+]=N%19)[C+]%19[C+2][C+](C%20=C([C+]=N%20)[C+]%20["
                "C+2][C+]%21C%22=C(N=[C+]%22)[C+]%22[C+2][C+]%23[C+2]C"
                "%24([C+2][C+](F)[C+]%22%24)C%22=C([C+]=N%22)C%22%24[C"
                "+2][C+]([C+2][C+](C%25=C(N=[C+]%25)C%25%26[C+2][C+](["
                "C+2][C+](C%27=C([C+]=N%27)C%27%28[C+2][C+]%29[C+2][C+"
                "](C%30=C(N=[C+]%30)[C+]%30[C+2][C+]([C+2]C%31([C+2][C"
                "+](F)[C+]%30%31)C%30=C([C+]=N%30)C%30%31[C+2][C+]([C+"
                "2][C+](C%32=C(N=[C+]%32)C%32%33[C+2][C+]([C+2][C+](C%"
                "34=C([C+]=N%34)C%34([C+2]%20)[C+2][C+](F)[C+]%21%34)["
                "C+]%32[C+](F)[C+2]%33)C%20=C7N=[C+]%20)[C+]%30[C+](F)"
                "[C+2]%31)C7=C(N=[C+]7)C7%20[C+2][C+]([C+2][C+](C%21=C"
                "%10[C+]=N%21)[C+]7[C+](F)[C+2]%20)C7=C(N=[C+]7)C7%10["
                "C+2][C+]([C+2][C+](C%20=C%14[C+]=N%20)[C+]7[C+](F)[C+"
                "2]%10)C7=C%29[C+]=N7)C7=C%23N=[C+]7)[C+]%27[C+](F)[C+"
                "2]%28)[C+]%25[C+](F)[C+2]%26)C7=C%17N=[C+]7)[C+]%22[C"
                "+](F)[C+2]%24)C7=C%18N=[C+]7)[C+]7[C+](F)[C+2]C7([C+2"
                "]%19)C7=C4[C+]=N7)[C+]2[C+](F)[C+2]%12)[C+]%15[C+](F)"
                "[C+2]%16)[C+]%11[C+](F)[C+2]%13)[C+]8[C+](F)[C+2]9)[C"
                "+]5[C+](F)[C+2]6)[C+]13"
            ),
            name=name,
        ),
    ),
)
def cage_twenty_plus_thirty(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
