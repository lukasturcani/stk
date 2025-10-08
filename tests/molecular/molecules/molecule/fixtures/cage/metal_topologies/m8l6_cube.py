import pytest

import stk

from ....case_data import CaseData
from ...building_blocks import get_iron_complex, get_tetratopic_linker


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.M8L6Cube(
                    building_blocks={
                        get_iron_complex(): range(8),
                        get_tetratopic_linker(): range(8, 14),
                    },
                ),
            ),
            smiles=(
                "[H]C1=C([H])C([H])=[N]2->[Fe+2]3456<-[N]7=C([H])C([H])=C([H]"
                ")C([H])=C7C([H])=[N]->3C3=C([H])C([H])=C(C7([H])C8=C([H])C(["
                "H])=C([N]9=C([H])C%10=C([H])C([H])=C([H])C([H])=[N]%10->[Fe+"
                "2]<-9%10%11%12<-[N]9=C([H])C([H])=C([H])C([H])=C9C([H])=[N]-"
                ">%10C9=C([H])C([H])=C(C([H])(C%10=C([H])C([H])=C([N]->4=C([H"
                "])C2=C1[H])C([H])=C%10[H])C1([H])C2=C([H])C([H])=C([N]4=C([H"
                "])C%10=C([H])C([H])=C([H])C([H])=[N]%10->[Fe+2]<-4%10%13%14<"
                "-[N]4=C([H])C([H])=C([H])C([H])=C4C([H])=[N]->%10C4=C([H])C("
                "[H])=C(C([H])(C%10=C([H])C([H])=C([N]->5=C([H])C5=C([H])C([H"
                "])=C([H])C([H])=[N]->65)C([H])=C%10[H])C5([H])C6=C([H])C([H]"
                ")=C([N]%10=C([H])C%15=C([H])C([H])=C([H])C([H])=[N]%15->[Fe+"
                "2]<-%10%15%16%17<-[N]%10=C([H])C([H])=C([H])C([H])=C%10C([H]"
                ")=[N]->%15C%10=C([H])C([H])=C(C7([H])C7=C([H])C([H])=C([N]%1"
                "5=C([H])C%18=C([H])C([H])=C([H])C([H])=[N]%18->[Fe+2]<-%15%1"
                "8%19%20<-[N]%15=C([H])C([H])=C([H])C([H])=C%15C([H])=[N]->%1"
                "8C%15=C([H])C([H])=C(C([H])(C%18=C([H])C([H])=C([N]->%16=C(["
                "H])C%16=C([H])C([H])=C([H])C([H])=[N]->%17%16)C([H])=C%18[H]"
                ")C%16([H])C%17=C([H])C([H])=C([N]%18=C([H])C%21=C([H])C([H])"
                "=C([H])C([H])=[N]%21->[Fe+2]<-%18%21%22(<-[N]%18=C([H])C([H]"
                ")=C([H])C([H])=C%18C([H])=[N]->%21C%18=C([H])C([H])=C5C([H])"
                "=C%18[H])<-[N]5=C([H])C([H])=C([H])C([H])=C5C([H])=[N]->%22C"
                "5=C([H])C([H])=C(C%18([H])C%21=C([H])C([H])=C([N]%22=C([H])C"
                "%23=C([H])C([H])=C([H])C([H])=[N]%23->[Fe+2]<-%22%23%24(<-[N"
                "]%22=C([H])C([H])=C([H])C([H])=C%22C([H])=[N]->%23C%22=C([H]"
                ")C([H])=C(C([H])(C%23=C([H])C([H])=C([N]->%19=C([H])C%19=C(["
                "H])C([H])=C([H])C([H])=[N]->%20%19)C([H])=C%23[H])C([H])(C%1"
                "9=C([H])C([H])=C([N]->%11=C([H])C%11=C([H])C([H])=C([H])C([H"
                "])=[N]->%12%11)C([H])=C%19[H])C%11=C([H])C([H])=C([N]%12=C(["
                "H])C%19=C([H])C([H])=C([H])C([H])=[N]%19->[Fe+2]<-%12%19%20("
                "<-[N]%12=C([H])C([H])=C([H])C([H])=C%12C([H])=[N]->%19C%12=C"
                "([H])C([H])=C1C([H])=C%12[H])<-[N]1=C([H])C([H])=C([H])C([H]"
                ")=C1C([H])=[N]->%20C1=C([H])C([H])=C(C%18([H])C%12=C([H])C(["
                "H])=C([N]->%13=C([H])C%13=C([H])C([H])=C([H])C([H])=[N]->%14"
                "%13)C([H])=C%12[H])C([H])=C1[H])C([H])=C%11[H])C([H])=C%22[H"
                "])<-[N]1=C([H])C([H])=C([H])C([H])=C1C([H])=[N]->%24C1=C([H]"
                ")C([H])=C%16C([H])=C1[H])C([H])=C%21[H])C([H])=C5[H])C([H])="
                "C%17[H])C([H])=C%15[H])C([H])=C7[H])C([H])=C%10[H])C([H])=C6"
                "[H])C([H])=C4[H])C([H])=C2[H])C([H])=C9[H])C([H])=C8[H])C([H"
                "])=C3[H]"
            ),
            name=name,
        ),
    ),
)
def metal_cage_m8l6_cube(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
