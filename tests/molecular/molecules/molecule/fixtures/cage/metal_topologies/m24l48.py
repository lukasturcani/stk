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
                "[H]C1=C([H])C2C3=C([H])C([H])=[N](->[Pd+2]45<-[N]6=C([H])C(["
                "H])=C(C7=C([H])C([H])=C([H])C(C8=C([H])C([H])=[N](->[Pd+2]9%1"
                "0<-[N]%11=C([H])C([H])=C(C%12=C([H])C([H])=C([H])C(C%13=C([H]"
                ")C([H])=[N](->[Pd+2]%14%15<-[N]%16=C([H])C([H])=C(C%17=C([H])"
                "C([H])=C([H])C(C%18=C([H])C([H])=[N](->[Pd+2]%19%20<-[N]%21=C"
                "([H])C([H])=C(C%22=C([H])C([H])=C([H])C(C%23=C([H])C([H])=[N]"
                "(->[Pd+2]%24%25<-[N]%26=C([H])C([H])=C(C%27=C([H])C([H])=C([H"
                "])C(C%28=C([H])C([H])=[N](->[Pd+2]%29%30<-[N]%31=C([H])C([H])"
                "=C(C%32=C([H])C([H])=C([H])C(C%33=C([H])C([H])=[N](->[Pd+2]"
                "%34%35<-[N]%36=C([H])C([H])=C(C%37=C([H])C(C%38=C([H])C([H])"
                "=[N](->[Pd+2]%39(<-[N]%40=C([H])C([H])=C(C%41=C([H])C([H])=C("
                "[H])C(C%42=C([H])C([H])=[N](->[Pd+2]%43%44<-[N]%45=C([H])C([H"
                "])=C(C%46=C([H])C([H])=C([H])C(C%47=C([H])C([H])=[N](->[Pd+2]"
                "%48%49<-[N]%50=C([H])C([H])=C(C%51=C([H])C([H])=C([H])C(C%52="
                "C([H])C([H])=[N](->[Pd+2](<-[N]%53=C([H])C([H])=C(C%54=C([H])"
                "C([H])=C([H])C(C%55=C([H])C([H])=[N]->%24C([H])=C%55[H])=C%54"
                "[H])C([H])=C%53[H])(<-[N]%24=C([H])C([H])=C(C%53=C([H])C([H]"
                ")=C([H])C(C%54=C([H])C([H])=[N](->[Pd+2]%55(<-[N]%56=C([H])C("
                "[H])=C(C%57=C([H])C(C%58=C([H])C([H])=[N]->%25C([H])=C%58[H])"
                "=C([H])C([H])=C%57[H])C([H])=C%56[H])<-[N]%25=C([H])C([H])=C("
                "C%56=C([H])C(C%57=C([H])C([H])=[N](->[Pd+2]%58(<-[N]%59=C([H]"
                ")C([H])=C(C%60=C([H])C(C%61=C([H])C([H])=[N]->%48C([H])=C%61["
                "H])=C([H])C([H])=C%60[H])C([H])=C%59[H])<-[N]%48=C([H])C([H])"
                "=C(C%59=C([H])C(C%60=C([H])C([H])=[N](->[Pd+2](<-[N]%61=C([H]"
                ")C([H])=C(C%62=C([H])C([H])=C([H])C(C%63=C([H])C([H])=[N](->["
                "Pd+2](<-[N]%64=C([H])C([H])=C(C%65=C([H])C([H])=C([H])C(C%66="
                "C([H])C([H])=[N]->%55C([H])=C%66[H])=C%65[H])C([H])=C%64[H])("
                "<-[N]%55=C([H])C([H])=C(C%64=C([H])C(C%65=C([H])C([H])=[N]"
                "->%14C([H])=C%65[H])=C([H])C([H])=C%64[H])C([H])=C%55[H])<-"
                "[N]%14=C([H])C([H])=C(C%55=C([H])C(C%64=C([H])C([H])=[N]"
                "->%19C([H])=C%64[H])=C([H])C([H])=C%55[H])C([H])=C%14[H])C("
                "[H])=C%63[H])=C%62[H])C([H])=C%61[H])(<-[N]%14=C([H])C([H])="
                "C(C%19=C([H])C(C%55=C([H])C([H])=[N]->9C([H])=C%55[H])=C([H]"
                ")C([H])=C%19[H])C([H])=C%14[H])<-[N]9=C([H])C([H])=C(C%14=C("
                "[H])C(C%19=C([H])C([H])=[N](->[Pd+2]%55(<-[N]%61=C([H])C([H]"
                ")=C(C%62=C([H])C([H])=C([H])C(C%63=C([H])C([H])=[N](->[Pd+2]"
                "(<-[N]%64=C([H])C([H])=C(C%65=C([H])C([H])=C([H])C(C%66=C([H]"
                ")C([H])=[N]->%58C([H])=C%66[H])=C%65[H])C([H])=C%64[H])(<-[N]"
                "%58=C([H])C([H])=C(C%64=C([H])C([H])=C([H])C(C%65=C([H])C([H]"
                ")=[N](->[Pd+2]%66(<-[N]%67=C([H])C([H])=C(C%68=C([H])C(C%69=C"
                "([H])C([H])=[N]->%43C([H])=C%69[H])=C([H])C([H])=C%68[H])C([H"
                "])=C%67[H])<-[N]%43=C([H])C([H])=C(C%67=C([H])C(C%68=C([H])C("
                "[H])=[N](->[Pd+2](<-[N]%69=C([H])C([H])=C(C%70=C([H])C([H])=C"
                "([H])C(C%71=C([H])C([H])=[N]->%55C([H])=C%71[H])=C%70[H])C([H"
                "])=C%69[H])(<-[N]%55=C([H])C([H])=C(C%69=C([H])C([H])=C([H])C"
                "(C%70=C([H])C([H])=[N](->[Pd+2](<-[N]%71=C([H])C([H])=C(C%72="
                "C([H])C([H])=C([H])C(C%73=C([H])C([H])=[N](->[Pd+2](<-[N]%74="
                "C([H])C([H])=C(C%75=C([H])C([H])=C([H])C(C%76=C([H])C([H])="
                "[N]->%66C([H])=C%76[H])=C%75[H])C([H])=C%74[H])(<-[N]%66=C("
                "[H])C([H])=C(C%74=C([H])C([H])=C([H])C(C%75=C([H])C([H])="
                "[N]->%34C([H])=C%75[H])=C%74[H])C([H])=C%66[H])<-[N]%34=C([H]"
                ")C([H])=C(C%66=C([H])C(C%74=C([H])C([H])=[N]->%44C([H])"
                "=C%74[H])=C([H])C([H])=C%66[H])C([H])=C%34[H])C([H])=C%73[H])"
                "=C%72[H])C([H])=C%71[H])(<-[N]%34=C([H])C([H])=C(C%44=C([H])"
                "C([H])=C([H])C(C%66=C([H])C([H])=[N](->[Pd+2]%71(<-[N]%72=C("
                "[H])C([H])=C(C%73=C([H])C([H])=C([H])C(C%74=C([H])C([H])="
                "[N]->%35C([H])=C%74[H])=C%73[H])C([H])=C%72[H])<-[N]%35=C("
                "[H])C([H])=C(C%72=C([H])C(C%73=C([H])C([H])=[N](->[Pd+2]"
                "(<-[N]%74=C([H])C([H])=C(C(=C1[H])C=2[H])C([H])=C%74[H])"
                "(<-[N]1=C([H])C([H])=C(C2=C([H])C([H])=C([H])C(C%74=C([H])"
                "C([H])=[N](->[Pd+2](<-[N]%75=C([H])C([H])=C(C%76=C([H])C("
                "[H])=C([H])C(C%77=C([H])C([H])=[N]->%20C([H])=C%77[H])"
                "=C%76[H])C([H])=C%75[H])(<-[N]%20=C([H])C([H])=C(C%75=C([H])"
                "C([H])=C([H])C(C%76=C([H])C([H])=[N]->%71C([H])=C%76[H])"
                "=C%75[H])C([H])=C%20[H])<-[N]%20=C([H])C([H])=C(C%71=C([H])"
                "C([H])=C([H])C(C%75=C([H])C([H])=[N]->%29C([H])=C%75[H])"
                "=C%71[H])C([H])=C%20[H])C([H])=C%74[H])=C2[H])C([H])=C1[H]"
                ")<-[N]1=C([H])C([H])=C(C2=C([H])C([H])=C([H])C(C%20=C([H])"
                "C([H])=[N]->%15C([H])=C%20[H])=C2[H])C([H])=C1[H])C([H])"
                "=C%73[H])=C([H])C([H])=C%72[H])C([H])=C%35[H])C([H])=C%66[H])"
                "=C%44[H])C([H])=C%34[H])<-[N]1=C([H])C([H])=C(C2=C([H])C"
                "(C%15=C([H])C([H])=[N]->4C([H])=C%15[H])=C([H])C([H])=C2[H])"
                "C([H])=C1[H])C([H])=C%70[H])=C%69[H])C([H])=C%55[H])<-[N]1=C"
                "([H])C([H])=C(C2=C([H])C(C4=C([H])C([H])=[N]->5C([H])=C4[H])"
                "=C([H])C([H])=C2[H])C([H])=C1[H])C([H])=C%68[H])=C([H])C([H])"
                "=C%67[H])C([H])=C%43[H])C([H])=C%65[H])=C%64[H])C([H])"
                "=C%58[H])<-[N]1=C([H])C([H])=C(C2=C([H])C(C4=C([H])C([H])="
                "[N]->%49C([H])=C4[H])=C([H])C([H])=C2[H])C([H])=C1[H])C([H]"
                ")=C%63[H])=C%62[H])C([H])=C%61[H])<-[N]1=C([H])C([H])=C(C2="
                "C([H])C(C4=C([H])C([H])=[N]->%10C([H])=C4[H])=C([H])C([H])="
                "C2[H])C([H])=C1[H])C([H])=C%19[H])=C([H])C([H])=C%14[H])C("
                "[H])=C9[H])C([H])=C%60[H])=C([H])C([H])=C%59[H])C([H])"
                "=C%48[H])C([H])=C%57[H])=C([H])C([H])=C%56[H])C([H])=C%25[H]"
                ")C([H])=C%54[H])=C%53[H])C([H])=C%24[H])<-[N]1=C([H])C([H])="
                "C(C2=C([H])C(C4=C([H])C([H])=[N]->%39C([H])=C4[H])=C([H])C"
                "([H])"
                "=C2[H])C([H])=C1[H])C([H])=C%52[H])=C%51[H])C([H])=C%50[H])C"
                "([H])=C%47[H])=C%46[H])C([H])=C%45[H])C([H])=C%42[H])=C%41[H]"
                ")C([H])=C%40[H])<-[N]1=C([H])C([H])=C(C2=C([H])C([H])=C([H])C"
                "(C4=C([H])C([H])=[N]->%30C([H])=C4[H])=C2[H])C([H])=C1[H])C("
                "[H])=C%38[H])=C([H])C([H])=C%37[H])C([H])=C%36[H])C([H])"
                "=C%33[H])=C%32[H])C([H])=C%31[H])C([H])=C%28[H])=C%27[H])C("
                "[H])=C%26[H])C([H])=C%23[H])=C%22[H])C([H])=C%21[H])C([H])="
                "C%18[H])=C%17[H])C([H])=C%16[H])C([H])=C%13[H])=C%12[H])C("
                "[H])=C%11[H])C([H])"
                "=C8[H])=C7[H])C([H])=C6[H])C([H])=C3[H]"
            ),
            name=name,
        ),
    ),
)
def metal_cage_m24l48(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
