import pytest
import stk

from ....case_data import CaseData
from ...building_blocks import get_linker, get_pd_atom


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                stk.cage.M12L24(
                    building_blocks=(
                        get_pd_atom(),
                        get_linker(),
                    ),
                    reaction_factory=stk.DativeReactionFactory(
                        stk.GenericReactionFactory(
                            bond_orders={
                                frozenset(
                                    {
                                        stk.GenericFunctionalGroup,
                                        stk.SingleAtom,
                                    }
                                ): 9
                            }
                        )
                    ),
                )
            ),
            smiles=(
                "[H]C1=C([H])C2=C([H])C(=C1[H])C1=C([H])C([H])=N(->["
                "Pd+2]34<-N5=C([H])C([H])=C(C([H])=C5[H])C5=C([H])C("
                "[H])=C([H])C(=C5[H])C5=C([H])C([H])=N(->[Pd+2]67<-N"
                "8=C([H])C([H])=C(C([H])=C8[H])C8=C([H])C([H])=C([H]"
                ")C(=C8[H])C8=C([H])C([H])=N(->[Pd+2]9(<-N%10=C([H])"
                "C([H])=C(C([H])=C%10[H])C%10=C([H])C([H])=C([H])C(="
                "C%10[H])C%10=C([H])C([H])=N(->[Pd+2]%11%12<-N%13=C("
                "[H])C([H])=C(C([H])=C%13[H])C%13=C([H])C(=C([H])C(["
                "H])=C%13[H])C%13=C([H])C([H])=N(->[Pd+2]%14(<-N%15="
                "C([H])C([H])=C(C([H])=C%15[H])C%15=C([H])C([H])=C(["
                "H])C(=C%15[H])C%15=C([H])C([H])=N(->[Pd+2]%16(<-N%1"
                "7=C([H])C([H])=C(C([H])=C%17[H])C%17=C([H])C([H])=C"
                "([H])C(=C%17[H])C%17=C([H])C([H])=N(->[Pd+2]%18(<-N"
                "%19=C([H])C([H])=C(C([H])=C%19[H])C%19=C([H])C(=C(["
                "H])C([H])=C%19[H])C%19=C([H])C([H])=N->%14C([H])=C%"
                "19[H])<-N%14=C([H])C([H])=C(C([H])=C%14[H])C%14=C(["
                "H])C(=C([H])C([H])=C%14[H])C%14=C([H])C([H])=N(->[P"
                "d+2]%19(<-N%20=C([H])C([H])=C(C([H])=C%20[H])C%20=C"
                "([H])C([H])=C([H])C(=C%20[H])C%20=C([H])C([H])=N(->"
                "[Pd+2]%21(<-N%22=C([H])C([H])=C(C([H])=C%22[H])C%22="
                "C([H])C([H])=C([H])C(=C%22[H])C%22=C([H])C([H])=N->%"
                "18C([H])=C%22[H])<-N%18=C([H])C([H])=C(C([H])=C%18[H"
                "])C%18=C([H])C(=C([H])C([H])=C%18[H])C%18=C([H])C([H"
                "])=N(->[Pd+2](<-N%22=C([H])C([H])=C(C([H])=C%22[H])C"
                "%22=C([H])C([H])=C([H])C(=C%22[H])C%22=C([H])C([H])="
                "N(->[Pd+2](<-N%23=C([H])C([H])=C(C([H])=C%23[H])C%2"
                "3=C([H])C([H])=C([H])C(=C%23[H])C%23=C([H])C([H])=N"
                "->%16C([H])=C%23[H])(<-N%16=C([H])C([H])=C(C([H])=C"
                "%16[H])C%16=C([H])C([H])=C([H])C(=C%16[H])C%16=C([H"
                "])C([H])=N->%21C([H])=C%16[H])<-N%16=C([H])C([H])=C"
                "2C([H])=C%16[H])C([H])=C%22[H])(<-N2=C([H])C([H])=C"
                "(C([H])=C2[H])C2=C([H])C([H])=C([H])C(=C2[H])C2=C(["
                "H])C([H])=N->6C([H])=C2[H])<-N2=C([H])C([H])=C(C([H"
                "])=C2[H])C2=C([H])C([H])=C([H])C(=C2[H])C2=C([H])C("
                "[H])=N(->[Pd+2](<-N6=C([H])C([H])=C(C([H])=C6[H])C6"
                "=C([H])C([H])=C([H])C(=C6[H])C6=C([H])C([H])=N->%11"
                "C([H])=C6[H])(<-N6=C([H])C([H])=C(C([H])=C6[H])C6=C"
                "([H])C(=C([H])C([H])=C6[H])C6=C([H])C([H])=N->%19C("
                "[H])=C6[H])<-N6=C([H])C([H])=C(C([H])=C6[H])C6=C([H]"
                ")C(=C([H])C([H])=C6[H])C6=C([H])C([H])=N->7C([H])=C6"
                "[H])C([H])=C2[H])C([H])=C%18[H])C([H])=C%20[H])<-N2="
                "C([H])C([H])=C(C([H])=C2[H])C2=C([H])C([H])=C([H])C("
                "=C2[H])C2=C([H])C([H])=N->%12C([H])=C2[H])C([H])=C%14"
                "[H])C([H])=C%17[H])<-N2=C([H])C([H])=C(C([H])=C2[H])"
                "C2=C([H])C(=C([H])C([H])=C2[H])C2=C([H])C([H])=N->3"
                "C([H])=C2[H])C([H])=C%15[H])<-N2=C([H])C([H])=C(C(["
                "H])=C2[H])C2=C([H])C([H])=C([H])C(=C2[H])C2=C([H])C"
                "([H])=N->9C([H])=C2[H])C([H])=C%13[H])C([H])=C%10[H"
                "])<-N2=C([H])C([H])=C(C([H])=C2[H])C2=C([H])C(=C([H"
                "])C([H])=C2[H])C2=C([H])C([H])=N->4C([H])=C2[H])C([H"
                "])=C8[H])C([H])=C5[H])C([H])=C1[H]"
            ),
            name=name,
        ),
    ),
)
def metal_cage_m12l24(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
