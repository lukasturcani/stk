import pytest

import stk

from ....case_data import CaseData
from ...building_blocks import (
    get_iron_complex,
    get_tetratopic_linker,
    get_tritopic_linker,
)


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.M6L2L3Prism(
                    building_blocks={
                        get_iron_complex(): range(6),
                        get_tritopic_linker(): range(6, 8),
                        get_tetratopic_linker(): range(8, 11),
                    },
                ),
            ),
            smiles=(
                "[H]C1=C([H])C([H])=[N]2->[Fe+2]3456<-[N]7=C([H])C([H])=C([H]"
                ")C([H])=C7C([H])=[N]->3C3=C([H])C([H])=C(C7([H])C8=C([H])C(["
                "H])=C([N]9=C([H])C%10=C([H])C([H])=C([H])C([H])=[N]%10->[Fe+"
                "2]<-9%10%11%12<-[N]9=C([H])C([H])=C([H])C([H])=C9C([H])=[N]-"
                ">%10C9=C([H])C([H])=C(N(C%10=C([H])C([H])=C([N]->4=C([H])C2="
                "C1[H])C([H])=C%10[H])C1=C([H])C([H])=C([N]2=C([H])C4=C([H])C"
                "([H])=C([H])C([H])=[N]4->[Fe+2]<-24%10%13<-[N]2=C([H])C([H])"
                "=C([H])C([H])=C2C([H])=[N]->4C2=C([H])C([H])=C(C([H])(C4=C(["
                "H])C([H])=C([N]->5=C([H])C5=C([H])C([H])=C([H])C([H])=[N]->6"
                "5)C([H])=C4[H])C4([H])C5=C([H])C([H])=C([N]6=C([H])C%14=C([H"
                "])C([H])=C([H])C([H])=[N]%14->[Fe+2]<-6%14%15%16<-[N]6=C([H]"
                ")C([H])=C([H])C([H])=C6C([H])=[N]->%14C6=C([H])C([H])=C(N%14"
                "C%17=C([H])C([H])=C([N]%18=C([H])C%19=C([H])C([H])=C([H])C(["
                "H])=[N]%19->[Fe+2]<-%18%19%20(<-[N]%18=C([H])C([H])=C([H])C("
                "[H])=C%18C([H])=[N]->%19C%18=C([H])C([H])=C(C7([H])C7=C([H])"
                "C([H])=C([N]->%15=C([H])C%15=C([H])C([H])=C([H])C([H])=[N]->"
                "%16%15)C([H])=C7[H])C([H])=C%18[H])<-[N]7=C([H])C([H])=C([H]"
                ")C([H])=C7C([H])=[N]->%20C7=C([H])C([H])=C(C([H])(C%15=C([H]"
                ")C([H])=C([N]->%11=C([H])C%11=C([H])C([H])=C([H])C([H])=[N]-"
                ">%12%11)C([H])=C%15[H])C([H])(C%11=C([H])C([H])=C([N]->%10=C"
                "([H])C%10=C([H])C([H])=C([H])C([H])=[N]->%13%10)C([H])=C%11["
                "H])C%10=C([H])C([H])=C([N]%11=C([H])C%12=C([H])C([H])=C([H])"
                "C([H])=[N]%12->[Fe+2]<-%11%12%13(<-[N]%11=C([H])C([H])=C([H]"
                ")C([H])=C%11C([H])=[N]->%12C%11=C([H])C([H])=C%14C([H])=C%11"
                "[H])<-[N]%11=C([H])C([H])=C([H])C([H])=C%11C([H])=[N]->%13C%"
                "11=C([H])C([H])=C4C([H])=C%11[H])C([H])=C%10[H])C([H])=C7[H]"
                ")C([H])=C%17[H])C([H])=C6[H])C([H])=C5[H])C([H])=C2[H])C([H]"
                ")=C1[H])C([H])=C9[H])C([H])=C8[H])C([H])=C3[H]"
            ),
            name=name,
        ),
    ),
)
def metal_cage_m6l2l3_prism(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
