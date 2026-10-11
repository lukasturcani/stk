import pytest

import stk

from ....case_data import CaseData
from ...building_blocks import get_linker, get_pd_atom


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.M9L18(
                    building_blocks=(get_pd_atom(), get_linker()),
                    reaction_factory=stk.DativeReactionFactory(
                        reaction_factory=stk.GenericReactionFactory(
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
                "[H]C1=C([H])C2C3=C([H])C([H])=[N](->[Pd+2]45<-[N]6=C([H])C([H"
                "])=C(C7=C([H])C([H])=C([H])C(C8=C([H])C([H])=[N](->[Pd+2]9%10"
                "<-[N]%11=C([H])C([H])=C(C%12=C([H])C([H])=C([H])C(C%13=C([H])"
                "C([H])=[N](->[Pd+2]%14%15<-[N]%16=C([H])C([H])=C(C%17=C([H])C"
                "([H])=C([H])C(C%18=C([H])C([H])=[N](->[Pd+2](<-[N]%19=C([H])C"
                "([H])=C(C%20=C([H])C(C%21=C([H])C([H])=[N]->4C([H])=C%21[H])="
                "C([H])C([H])=C%20[H])C([H])=C%19[H])(<-[N]4=C([H])C([H])=C(C%"
                "19=C([H])C(C%20=C([H])C([H])=[N]->9C([H])=C%20[H])=C([H])C([H"
                "])=C%19[H])C([H])=C4[H])<-[N]4=C([H])C([H])=C(C9=C([H])C(C%19"
                "=C([H])C([H])=[N](->[Pd+2]%20(<-[N]%21=C([H])C([H])=C(C%22=C("
                "[H])C(C%23=C([H])C([H])=[N](->[Pd+2]%24(<-[N]%25=C([H])C([H])"
                "=C(C%26=C([H])C([H])=C([H])C(C%27=C([H])C([H])=[N](->[Pd+2]%2"
                "8(<-[N]%29=C([H])C([H])=C(C%30=C([H])C([H])=C([H])C(C%31=C([H"
                "])C([H])=[N](->[Pd+2](<-[N]%32=C([H])C([H])=C(C%33=C([H])C([H"
                "])=C([H])C(C%34=C([H])C([H])=[N]->%20C([H])=C%34[H])=C%33[H])"
                "C([H])=C%32[H])(<-[N]%20=C([H])C([H])=C(C%32=C([H])C(C%33=C(["
                "H])C([H])=[N](->[Pd+2](<-[N]%34=C([H])C([H])=C(C(=C1[H])C=2[H"
                "])C([H])=C%34[H])(<-[N]1=C([H])C([H])=C(C2=C([H])C([H])=C([H]"
                ")C(C%34=C([H])C([H])=[N]->%28C([H])=C%34[H])=C2[H])C([H])=C1["
                "H])<-[N]1=C([H])C([H])=C(C2=C([H])C([H])=C([H])C(C%28=C([H])C"
                "([H])=[N]->%10C([H])=C%28[H])=C2[H])C([H])=C1[H])C([H])=C%33["
                "H])=C([H])C([H])=C%32[H])C([H])=C%20[H])<-[N]1=C([H])C([H])=C"
                "(C2=C([H])C(C%10=C([H])C([H])=[N]->%24C([H])=C%10[H])=C([H])C"
                "([H])=C2[H])C([H])=C1[H])C([H])=C%31[H])=C%30[H])C([H])=C%29["
                "H])<-[N]1=C([H])C([H])=C(C2=C([H])C([H])=C([H])C(C%10=C([H])C"
                "([H])=[N]->%14C([H])=C%10[H])=C2[H])C([H])=C1[H])C([H])=C%27["
                "H])=C%26[H])C([H])=C%25[H])<-[N]1=C([H])C([H])=C(C2=C([H])C(["
                "H])=C([H])C(C%10=C([H])C([H])=[N]->%15C([H])=C%10[H])=C2[H])C"
                "([H])=C1[H])C([H])=C%23[H])=C([H])C([H])=C%22[H])C([H])=C%21["
                "H])<-[N]1=C([H])C([H])=C(C2=C([H])C(C%10=C([H])C([H])=[N]->5C"
                "([H])=C%10[H])=C([H])C([H])=C2[H])C([H])=C1[H])C([H])=C%19[H]"
                ")=C([H])C([H])=C9[H])C([H])=C4[H])C([H])=C%18[H])=C%17[H])C(["
                "H])=C%16[H])C([H])=C%13[H])=C%12[H])C([H])=C%11[H])C([H])=C8["
                "H])=C7[H])C([H])=C6[H])C([H])=C3[H]"
            ),
            name=name,
        ),
    ),
)
def metal_cage_m9l18(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
