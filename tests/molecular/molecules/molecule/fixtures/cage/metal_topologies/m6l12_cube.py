import pytest

import stk

from ....case_data import CaseData
from ...building_blocks import get_linker, get_pd_atom


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.M6L12Cube(
                    building_blocks={
                        get_pd_atom(): range(6),
                        get_linker(): range(6, 18),
                    },
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
                "[H]C1=C([H])C2C3=C([H])C([H])=[N](->[Pd+2]45<-[N]6=C([H])C(["
                "H])=C(C7=C([H])C([H])=C([H])C(C8=C([H])C([H])=[N](->[Pd+2]9%"
                "10<-[N]%11=C([H])C([H])=C(C%12=C([H])C([H])=C([H])C(C%13=C(["
                "H])C([H])=[N](->[Pd+2]%14%15<-[N]%16=C([H])C([H])=C(C%17=C(["
                "H])C([H])=C([H])C(C%18=C([H])C([H])=[N](->[Pd+2](<-[N]%19=C("
                "[H])C([H])=C(C%20=C([H])C(C%21=C([H])C([H])=[N](->[Pd+2](<-["
                "N]%22=C([H])C([H])=C(C(=C1[H])C=2[H])C([H])=C%22[H])(<-[N]1="
                "C([H])C([H])=C(C2=C([H])C([H])=C([H])C(C%22=C([H])C([H])=[N]"
                "->%14C([H])=C%22[H])=C2[H])C([H])=C1[H])<-[N]1=C([H])C([H])="
                "C(C2=C([H])C([H])=C([H])C(C%14=C([H])C([H])=[N](->[Pd+2](<-["
                "N]%22=C([H])C([H])=C(C%23=C([H])C(C%24=C([H])C([H])=[N]->4C("
                "[H])=C%24[H])=C([H])C([H])=C%23[H])C([H])=C%22[H])(<-[N]4=C("
                "[H])C([H])=C(C%22=C([H])C(C%23=C([H])C([H])=[N]->9C([H])=C%2"
                "3[H])=C([H])C([H])=C%22[H])C([H])=C4[H])<-[N]4=C([H])C([H])="
                "C(C9=C([H])C(C%22=C([H])C([H])=[N]->%15C([H])=C%22[H])=C([H]"
                ")C([H])=C9[H])C([H])=C4[H])C([H])=C%14[H])=C2[H])C([H])=C1[H"
                "])C([H])=C%21[H])=C([H])C([H])=C%20[H])C([H])=C%19[H])(<-[N]"
                "1=C([H])C([H])=C(C2=C([H])C(C4=C([H])C([H])=[N]->5C([H])=C4["
                "H])=C([H])C([H])=C2[H])C([H])=C1[H])<-[N]1=C([H])C([H])=C(C2"
                "=C([H])C(C4=C([H])C([H])=[N]->%10C([H])=C4[H])=C([H])C([H])="
                "C2[H])C([H])=C1[H])C([H])=C%18[H])=C%17[H])C([H])=C%16[H])C("
                "[H])=C%13[H])=C%12[H])C([H])=C%11[H])C([H])=C8[H])=C7[H])C(["
                "H])=C6[H])C([H])=C3[H]"
            ),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.M6L12Cube(
                    building_blocks={
                        get_pd_atom(): range(6),
                        get_linker(): range(6, 18),
                    },
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
                    scale_multiplier=1.2,
                ),
            ),
            smiles=(
                "[H]C1=C([H])C2C3=C([H])C([H])=[N](->[Pd+2]45<-[N]6=C([H])C(["
                "H])=C(C7=C([H])C([H])=C([H])C(C8=C([H])C([H])=[N](->[Pd+2]9%"
                "10<-[N]%11=C([H])C([H])=C(C%12=C([H])C([H])=C([H])C(C%13=C(["
                "H])C([H])=[N](->[Pd+2]%14%15<-[N]%16=C([H])C([H])=C(C%17=C(["
                "H])C([H])=C([H])C(C%18=C([H])C([H])=[N](->[Pd+2](<-[N]%19=C("
                "[H])C([H])=C(C%20=C([H])C(C%21=C([H])C([H])=[N](->[Pd+2](<-["
                "N]%22=C([H])C([H])=C(C(=C1[H])C=2[H])C([H])=C%22[H])(<-[N]1="
                "C([H])C([H])=C(C2=C([H])C([H])=C([H])C(C%22=C([H])C([H])=[N]"
                "->%14C([H])=C%22[H])=C2[H])C([H])=C1[H])<-[N]1=C([H])C([H])="
                "C(C2=C([H])C([H])=C([H])C(C%14=C([H])C([H])=[N](->[Pd+2](<-["
                "N]%22=C([H])C([H])=C(C%23=C([H])C(C%24=C([H])C([H])=[N]->4C("
                "[H])=C%24[H])=C([H])C([H])=C%23[H])C([H])=C%22[H])(<-[N]4=C("
                "[H])C([H])=C(C%22=C([H])C(C%23=C([H])C([H])=[N]->9C([H])=C%2"
                "3[H])=C([H])C([H])=C%22[H])C([H])=C4[H])<-[N]4=C([H])C([H])="
                "C(C9=C([H])C(C%22=C([H])C([H])=[N]->%15C([H])=C%22[H])=C([H]"
                ")C([H])=C9[H])C([H])=C4[H])C([H])=C%14[H])=C2[H])C([H])=C1[H"
                "])C([H])=C%21[H])=C([H])C([H])=C%20[H])C([H])=C%19[H])(<-[N]"
                "1=C([H])C([H])=C(C2=C([H])C(C4=C([H])C([H])=[N]->5C([H])=C4["
                "H])=C([H])C([H])=C2[H])C([H])=C1[H])<-[N]1=C([H])C([H])=C(C2"
                "=C([H])C(C4=C([H])C([H])=[N]->%10C([H])=C4[H])=C([H])C([H])="
                "C2[H])C([H])=C1[H])C([H])=C%18[H])=C%17[H])C([H])=C%16[H])C("
                "[H])=C%13[H])=C%12[H])C([H])=C%11[H])C([H])=C8[H])=C7[H])C(["
                "H])=C6[H])C([H])=C3[H]"
            ),
            name=name,
        ),
    ),
)
def metal_cage_m6l12_cube(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
