import pytest

import stk

from ....case_data import CaseData
from ...building_blocks import get_linker, get_pd_atom


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.M24L48(
                    building_blocks={
                        get_pd_atom(): range(24),
                        get_linker(): range(24, 72),
                    },
                    reaction_factory=stk.DativeReactionFactory(
                        stk.GenericReactionFactory(
                            bond_orders={
                                frozenset(
                                    {
                                        stk.GenericFunctionalGroup,
                                        stk.SingleAtom,
                                    }
                                ): 9,
                            },
                        ),
                    ),
                ),
            ),
            smiles=(
                "[H]C1=C([H])C2C3=C([H])C([H])=[N](->[Pd+2]45<-[N]6=C([H])C("
                "[H])=C(C7=C([H])C([H])=C([H])C(C8=C([H])C([H])=[N](->[Pd+2]"
                "9%10<-[N]%11=C([H])C([H])=C(C%12=C([H])C([H])=C([H])C(C%13="
                "C([H])C([H])=[N](->[Pd+2]%14%15<-[N]%16=C([H])C([H])=C(C%17"
                "=C([H])C([H])=C([H])C(C%18=C([H])C([H])=[N](->[Pd+2]%19%20<"
                "-[N]%21=C([H])C([H])=C(C%22=C([H])C([H])=C([H])C(C%23=C([H]"
                ")C([H])=[N](->[Pd+2]%24%25<-[N]%26=C([H])C([H])=C(C%27=C([H"
                "])C([H])=C([H])C(C%28=C([H])C([H])=[N](->[Pd+2]%29%30<-[N]%"
                "31=C([H])C([H])=C(C%32=C([H])C(C%33=C([H])C([H])=[N](->[Pd+"
                "2]%34%35<-[N]%36=C([H])C([H])=C(C%37=C([H])C([H])=C([H])C(C"
                "%38=C([H])C([H])=[N](->[Pd+2](<-[N]%39=C([H])C([H])=C(C%40="
                "C([H])C([H])=C([H])C(C%41=C([H])C([H])=[N]->%19C([H])=C%41["
                "H])=C%40[H])C([H])=C%39[H])(<-[N]%19=C([H])C([H])=C(C%39=C("
                "[H])C([H])=C([H])C(C%40=C([H])C([H])=[N]->%24C([H])=C%40[H]"
                ")=C%39[H])C([H])=C%19[H])<-[N]%19=C([H])C([H])=C(C%24=C([H]"
                ")C(C%39=C([H])C([H])=[N](->[Pd+2]%40%41<-[N]%42=C([H])C([H]"
                ")=C(C%43=C([H])C([H])=C([H])C(C%44=C([H])C([H])=[N](->[Pd+2"
                "]%45%46<-[N]%47=C([H])C([H])=C(C%48=C([H])C([H])=C([H])C(C%"
                "49=C([H])C([H])=[N](->[Pd+2]%50%51<-[N]%52=C([H])C([H])=C(C"
                "%53=C([H])C(C%54=C([H])C([H])=[N](->[Pd+2](<-[N]%55=C([H])C"
                "([H])=C(C%56=C([H])C([H])=C([H])C(C%57=C([H])C([H])=[N](->["
                "Pd+2]%58(<-[N]%59=C([H])C([H])=C(C%60=C([H])C([H])=C([H])C("
                "C%61=C([H])C([H])=[N](->[Pd+2]%62(<-[N]%63=C([H])C([H])=C(C"
                "%64=C([H])C([H])=C([H])C(C%65=C([H])C([H])=[N](->[Pd+2]%66("
                "<-[N]%67=C([H])C([H])=C(C%68=C([H])C([H])=C([H])C(C%69=C([H"
                "])C([H])=[N](->[Pd+2]%70(<-[N]%71=C([H])C([H])=C(C%72=C([H]"
                ")C([H])=C([H])C(C%73=C([H])C([H])=[N]->%45C([H])=C%73[H])=C"
                "%72[H])C([H])=C%71[H])<-[N]%45=C([H])C([H])=C(C%71=C([H])C("
                "C%72=C([H])C([H])=[N](->[Pd+2](<-[N]%73=C([H])C([H])=C(C%74"
                "=C([H])C([H])=C([H])C(C%75=C([H])C([H])=[N](->[Pd+2](<-[N]%"
                "76=C([H])C([H])=C(C%77=C([H])C([H])=C([H])C(C%78=C([H])C([H"
                "])=[N]->%14C([H])="
                "C%78[H])=C%77[H])C([H])=C%76[H])(<-[N]%14=C([H])C([H])=C(C%"
                "76=C([H])C([H])=C([H])C(C%77=C([H])C([H])=[N]->%70C([H])=C%"
                "77[H])=C%76[H])C([H])=C%14[H])<-[N]%14=C([H])C([H])=C(C%70="
                "C([H])C(C%76=C([H])C([H])=[N](->[Pd+2](<-[N]%77=C([H])C([H]"
                ")=C(C%78=C([H])C([H])=C([H])C(C%79=C([H])C([H])=[N]->%40C(["
                "H])=C%79[H])=C%78[H])C([H])=C%77[H])(<-[N]%40=C([H])C([H])="
                "C(C%77=C([H])C([H])=C([H])C(C%78=C([H])C([H])=[N]->%20C([H]"
                ")=C%78[H])=C%77[H])C([H])=C%40[H])<-[N]%20=C([H])C([H])=C(C"
                "%40=C([H])C([H])=C([H])C"
                "(C%77=C([H])C([H])=[N]->%46C([H])=C%77[H])=C%40[H])C([H])=C"
                "%20[H])C([H])=C%76[H])=C([H])C([H])=C%70[H])C([H])=C%14[H])"
                "C([H])=C%75[H])=C%74[H])C([H])=C%73[H])(<-[N]%14=C([H])C([H"
                "])=C(C%20=C([H])C([H])=C([H])C(C%40=C([H])C([H])=[N]->9C([H"
                "])=C%40[H])=C%20[H])C([H])=C%14[H])<-[N]9=C([H])C([H])=C(C%"
                "14=C([H])C(C%20=C([H])C([H])=[N](->[Pd+2](<-[N]%40=C([H])C("
                "[H])=C(C(=C1[H])C=2[H])C([H])=C%40[H])(<-[N]1=C([H])C([H])="
                "C(C2=C([H])C([H])=C([H])C(C%40=C([H])C([H])=[N]->%62C([H])="
                "C%40[H])=C2[H])C([H])=C1"
                "[H])<-[N]1=C([H])C([H])=C(C2=C([H])C([H])=C([H])C(C%40=C([H"
                "])C([H])=[N]->%66C([H])=C%40[H])=C2[H])C([H])=C1[H])C([H])="
                "C%20[H])=C([H])C([H])=C%14[H])C([H])=C9[H])C([H])=C%72[H])="
                "C([H])C([H])=C%71[H])C([H])=C%45[H])C([H])=C%69[H])=C%68[H]"
                ")C([H])=C%67[H])<-[N]1=C([H])C([H])=C(C2=C([H])C([H])=C([H]"
                ")C(C9=C([H])C([H])=[N]->%50C([H])=C9[H])=C2[H])C([H])=C1[H]"
                ")C([H])=C%65[H])=C%64[H])C([H])=C%63[H])<-[N]1=C([H])C([H])"
                "=C(C2=C([H])C([H])=C([H])C(C9=C([H])C([H])=[N](->[Pd+2](<-["
                "N]%14=C([H])C([H])=C(C%2"
                "0=C([H])C([H])=C([H])C(C%40=C([H])C([H])=[N](->[Pd+2](<-[N]"
                "%45=C([H])C([H])=C(C%46=C([H])C([H])=C([H])C(C%50=C([H])C(["
                "H])=[N](->[Pd+2](<-[N]%62=C([H])C([H])=C(C%63=C([H])C([H])="
                "C([H])C(C%64=C([H])C([H])=[N]->%25C([H])=C%64[H])=C%63[H])C"
                "([H])=C%62[H])(<-[N]%25=C([H])C([H])=C(C%62=C([H])C(C%63=C("
                "[H])C([H])=[N]->%10C([H])=C%63[H])=C([H])C([H])=C%62[H])C(["
                "H])=C%25[H])<-[N]%10=C([H])C([H])=C(C%25=C([H])C(C%62=C([H]"
                ")C([H])=[N]->%15C([H])=C%62[H])=C([H])C([H])=C%25[H])C([H])"
                "=C%10[H])C([H])=C%50[H"
                "])=C%46[H])C([H])=C%45[H])(<-[N]%10=C([H])C([H])=C(C%15=C(["
                "H])C([H])=C([H])C(C%25=C([H])C([H])=[N]->%29C([H])=C%25[H])"
                "=C%15[H])C([H])=C%10[H])<-[N]%10=C([H])C([H])=C(C%15=C([H])"
                "C(C%25=C([H])C([H])=[N]->4C([H])=C%25[H])=C([H])C([H])=C%15"
                "[H])C([H])=C%10[H])C([H])=C%40[H])=C%20[H])C([H])=C%14[H])("
                "<-[N]4=C([H])C([H])=C(C%10=C([H])C([H])=C([H])C(C%14=C([H])"
                "C([H])=[N](->[Pd+2](<-[N]%15=C([H])C([H])=C(C%20=C([H])C([H"
                "])=C([H])C(C%25=C([H])C([H])=[N]->%30C([H])=C%25[H])=C%20[H"
                "])C([H])=C%15[H])(<-[N]%"
                "15=C([H])C([H])=C(C%20=C([H])C(C%25=C([H])C([H])=[N]->%58C("
                "[H])=C%25[H])=C([H])C([H])=C%20[H])C([H])=C%15[H])<-[N]%15="
                "C([H])C([H])=C(C%20=C([H])C(C%25=C([H])C([H])=[N]->%34C([H]"
                ")=C%25[H])=C([H])C([H])=C%20[H])C([H])=C%15[H])C([H])=C%14["
                "H])=C%10[H])C([H])=C4[H])<-[N]4=C([H])C([H])=C(C%10=C([H])C"
                "(C%14=C([H])C([H])=[N]->5C([H])=C%14[H])=C([H])C([H])=C%10["
                "H])C([H])=C4[H])C([H])=C9[H])=C2[H])C([H])=C1[H])C([H])=C%6"
                "1[H])=C%60[H])C([H])=C%59[H])<-[N]1=C([H])C([H])=C(C2=C([H]"
                ")C([H])=C([H])C(C4=C"
                "([H])C([H])=[N]->%51C([H])=C4[H])=C2[H])C([H])=C1[H])C([H])"
                "=C%57[H])=C%56[H])C([H])=C%55[H])(<-[N]1=C([H])C([H])=C(C2="
                "C([H])C([H])=C([H])C(C4=C([H])C([H])=[N]->%35C([H])=C4[H])="
                "C2[H])C([H])=C1[H])<-[N]1=C([H])C([H])=C(C2=C([H])C([H])=C("
                "[H])C(C4=C([H])C([H])=[N]->%41C([H])=C4[H])=C2[H])C([H])=C1"
                "[H])C([H])=C%54[H])=C([H])C([H])=C%53[H])C([H])=C%52[H])C(["
                "H])=C%49[H])=C%48[H])C([H])=C%47[H])C([H])=C%44[H])=C%43[H]"
                ")C([H])=C%42[H])C([H])=C%39[H])=C([H])C([H])=C%24[H])C([H])"
                "=C%19[H])C([H])=C%38[H"
                "])=C%37[H])C([H])=C%36[H])C([H])=C%33[H])=C([H])C([H])=C%32"
                "[H])C([H])=C%31[H])C([H])=C%28[H])=C%27[H])C([H])=C%26[H])C"
                "([H])=C%23[H])=C%22[H])C([H])=C%21[H])C([H])=C%18[H])=C%17["
                "H])C([H])=C%16[H])C([H])=C%13[H])=C%12[H])C([H])=C%11[H])C("
                "[H])=C8[H])=C7[H])C([H])=C6[H])C([H])=C3[H]"
            ),
            name=name,
        ),
    ),
)
def metal_cage_m24l48(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
