import pytest
import stk

from ....case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.EightPlusTwelve(
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
                "2][C+]7[C+2][C+](C8=C(N=[C+]8)C89[C+2][C+]%10[C+2][C+"
                "](C%11=C([C+]=N%11)C%11%12[C+2][C+]([C+2][C+](C%13=C2"
                "[C+]=N%13)[C+]%11[C+](F)[C+2]%12)C2=C(N=[C+]2)[C+]2[C"
                "+2][C+]%11C%12=C(N=[C+]%12)C%12%13[C+2][C+]([C+2][C+]"
                "(C%14=C([C+]=N%14)[C+]%14[C+2][C+](C%15=C7[C+]=N%15)["
                "C+]7[C+](F)[C+2]C7([C+2]%14)C7=C([C+]=N7)[C+]7[C+2][C"
                "+]([C+2]C%14([C+2][C+](F)[C+]7%14)C7=C4[C+]=N7)C4=C(N"
                "=[C+]4)C4([C+2]2)[C+2][C+](F)[C+]%114)[C+]%12[C+](F)["
                "C+2]%13)C2=C%10[C+]=N2)[C+]8[C+](F)[C+2]9)[C+]5[C+](F"
                ")[C+2]6)[C+]13"
            ),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.EightPlusTwelve(
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
                    optimizer=stk.MCHammer(
                        num_steps=150,
                        random_seed=1000,
                    ),
                ),
            ),
            smiles=(
                "F[C+]1[C+2]C23[C+2][C+]4[C+2][C+](C5=C(N=[C+]5)C56[C+"
                "2][C+]7[C+2][C+](C8=C(N=[C+]8)C89[C+2][C+]%10[C+2][C+"
                "](C%11=C([C+]=N%11)C%11%12[C+2][C+]([C+2][C+](C%13=C2"
                "[C+]=N%13)[C+]%11[C+](F)[C+2]%12)C2=C(N=[C+]2)[C+]2[C"
                "+2][C+]%11C%12=C(N=[C+]%12)C%12%13[C+2][C+]([C+2][C+]"
                "(C%14=C([C+]=N%14)[C+]%14[C+2][C+](C%15=C7[C+]=N%15)["
                "C+]7[C+](F)[C+2]C7([C+2]%14)C7=C([C+]=N7)[C+]7[C+2][C"
                "+]([C+2]C%14([C+2][C+](F)[C+]7%14)C7=C4[C+]=N7)C4=C(N"
                "=[C+]4)C4([C+2]2)[C+2][C+](F)[C+]%114)[C+]%12[C+](F)["
                "C+2]%13)C2=C%10[C+]=N2)[C+]8[C+](F)[C+2]9)[C+]5[C+](F"
                ")[C+2]6)[C+]13"
            ),
            name=name,
        ),
    ),
)
def cage_eight_plus_twelve(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
