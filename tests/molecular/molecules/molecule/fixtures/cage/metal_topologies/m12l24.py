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
                "[H]C1=C([H])C2C3=C([H])C([H])=[N](->[Pd+2]45<-[N]6=C([H])C([H"
                "])=C(C7=C([H])C([H])=C([H])C(C8=C([H])C([H])=[N](->[Pd+2]9%10"
                "<-[N]%11=C([H])C([H])=C(C%12=C([H])C([H])=C([H])C(C%13=C([H])"
                "C([H])=[N](->[Pd+2]%14(<-[N]%15=C([H])C([H])=C(C%16=C([H])C(C"
                "%17=C([H])C([H])=[N]->4C([H])=C%17[H])=C([H])C([H])=C%16[H])C"
                "([H])=C%15[H])<-[N]4=C([H])C([H])=C(C%15=C([H])C(C%16=C([H])C"
                "([H])=[N](->[Pd+2]%17%18<-[N]%19=C([H])C([H])=C(C%20=C([H])C("
                "[H])=C([H])C(C%21=C([H])C([H])=[N](->[Pd+2]%22(<-[N]%23=C([H]"
                ")C([H])=C(C%24=C([H])C([H])=C([H])C(C%25=C([H])C([H])=[N]->%1"
                "4C([H])=C%25[H])=C%24[H])C([H])=C%23[H])<-[N]%14=C([H])C([H])"
                "=C(C%23=C([H])C([H])=C([H])C(C%24=C([H])C([H])=[N](->[Pd+2]%2"
                "5(<-[N]%26=C([H])C([H])=C(C%27=C([H])C(C%28=C([H])C([H])=[N]("
                "->[Pd+2]%29(<-[N]%30=C([H])C([H])=C(C%31=C([H])C([H])=C([H])C"
                "(C%32=C([H])C([H])=[N]->9C([H])=C%32[H])=C%31[H])C([H])=C%30["
                "H])<-[N]9=C([H])C([H])=C(C%30=C([H])C([H])=C([H])C(C%31=C([H]"
                ")C([H])=[N](->[Pd+2]%32%33<-[N]%34=C([H])C([H])=C(C%35=C([H])"
                "C(C%36=C([H])C([H])=[N](->[Pd+2](<-[N]%37=C([H])C([H])=C(C(=C"
                "1[H])C=2[H])C([H])=C%37[H])(<-[N]1=C([H])C([H])=C(C2=C([H])C("
                "[H])=C([H])C(C%37=C([H])C([H])=[N]->%29C([H])=C%37[H])=C2[H])"
                "C([H])=C1[H])<-[N]1=C([H])C([H])=C(C2=C([H])C([H])=C([H])C(C%"
                "29=C([H])C([H])=[N](->[Pd+2](<-[N]%37=C([H])C([H])=C(C%38=C(["
                "H])C(C%39=C([H])C([H])=[N](->[Pd+2](<-[N]%40=C([H])C([H])=C(C"
                "%41=C([H])C([H])=C([H])C(C%42=C([H])C([H])=[N]->%17C([H])=C%4"
                "2[H])=C%41[H])C([H])=C%40[H])(<-[N]%17=C([H])C([H])=C(C%40=C("
                "[H])C([H])=C([H])C(C%41=C([H])C([H])=[N](->[Pd+2](<-[N]%42=C("
                "[H])C([H])=C(C%43=C([H])C([H])=C([H])C(C%44=C([H])C([H])=[N]-"
                ">%22C([H])=C%44[H])=C%43[H])C([H])=C%42[H])(<-[N]%22=C([H])C("
                "[H])=C(C%42=C([H])C([H])=C([H])C(C%43=C([H])C([H])=[N]->%32C("
                "[H])=C%43[H])=C%42[H])C([H])=C%22[H])<-[N]%22=C([H])C([H])=C("
                "C%32=C([H])C([H])=C([H])C(C%42=C([H])C([H])=[N]->%25C([H])=C%"
                "42[H])=C%32[H])C([H])=C%22[H])C([H])=C%41[H])=C%40[H])C([H])="
                "C%17[H])<-[N]%17=C([H])C([H])=C(C%22=C([H])C([H])=C([H])C(C%2"
                "5=C([H])C([H])=[N]->%33C([H])=C%25[H])=C%22[H])C([H])=C%17[H]"
                ")C([H])=C%39[H])=C([H])C([H])=C%38[H])C([H])=C%37[H])(<-[N]%1"
                "7=C([H])C([H])=C(C%22=C([H])C(C%25=C([H])C([H])=[N]->5C([H])="
                "C%25[H])=C([H])C([H])=C%22[H])C([H])=C%17[H])<-[N]5=C([H])C(["
                "H])=C(C%17=C([H])C(C%22=C([H])C([H])=[N]->%18C([H])=C%22[H])="
                "C([H])C([H])=C%17[H])C([H])=C5[H])C([H])=C%29[H])=C2[H])C([H]"
                ")=C1[H])C([H])=C%36[H])=C([H])C([H])=C%35[H])C([H])=C%34[H])C"
                "([H])=C%31[H])=C%30[H])C([H])=C9[H])C([H])=C%28[H])=C([H])C(["
                "H])=C%27[H])C([H])=C%26[H])<-[N]1=C([H])C([H])=C(C2=C([H])C(C"
                "5=C([H])C([H])=[N]->%10C([H])=C5[H])=C([H])C([H])=C2[H])C([H]"
                ")=C1[H])C([H])=C%24[H])=C%23[H])C([H])=C%14[H])C([H])=C%21[H]"
                ")=C%20[H])C([H])=C%19[H])C([H])=C%16[H])=C([H])C([H])=C%15[H]"
                ")C([H])=C4[H])C([H])=C%13[H])=C%12[H])C([H])=C%11[H])C([H])=C"
                "8[H])=C7[H])C([H])=C6[H])C([H])=C3[H]"
            ),
            name=name,
        ),
    ),
)
def metal_cage_m12l24(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
