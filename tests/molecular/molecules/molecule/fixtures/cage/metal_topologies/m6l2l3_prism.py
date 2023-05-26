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
                "[H]C1=C([H])C([H])=N2->[Fe+2]3456<-N7=C([H])C([H])=C("
                "[H])C([H])=C7C([H])=N->3C3=C([H])C([H])=C(C([H])=C3[H"
                "])C3([H])C7=C([H])C([H])=C(C([H])=C7[H])N7->[Fe+2]89%"
                "10(<-N%11=C([H])C([H])=C([H])C([H])=C%11C([H])=N->8C8"
                "=C([H])C([H])=C(C([H])=C8[H])N(C8=C([H])C([H])=C(C([H"
                "])=C8[H])N->4=C([H])C2=C1[H])C1=C([H])C([H])=C(C([H])"
                "=C1[H])N1->[Fe+2]248(<-N%11=C([H])C([H])=C([H])C([H])"
                "=C%11C=1[H])<-N1=C([H])C([H])=C([H])C([H])=C1C([H])=N"
                "->2C1=C([H])C([H])=C(C([H])=C1[H])C([H])(C1=C([H])C(["
                "H])=C(C([H])=C1[H])N->5=C([H])C1=C([H])C([H])=C([H])C"
                "([H])=N->61)C1([H])C2=C([H])C([H])=C(C([H])=C2[H])N2-"
                ">[Fe+2]56%11(<-N%12=C([H])C([H])=C([H])C([H])=C%12C(["
                "H])=N->5C5=C([H])C([H])=C(C([H])=C5[H])N5C%12=C([H])C"
                "([H])=C(C([H])=C%12[H])N%12->[Fe+2]%13%14(<-N%15=C([H"
                "])C([H])=C([H])C([H])=C%15C=%12[H])(<-N%12=C([H])C([H"
                "])=C([H])C([H])=C%12C([H])=N->%13C%12=C([H])C([H])=C("
                "C([H])=C%12[H])C3([H])C3=C([H])C([H])=C(C([H])=C3[H])"
                "N3->[Fe+2]%12%13(<-N%15=C([H])C([H])=C([H])C([H])=C%1"
                "5C([H])=N->%12C%12=C([H])C([H])=C5C([H])=C%12[H])(<-N"
                "5=C([H])C([H])=C([H])C([H])=C5C=3[H])<-N3=C([H])C([H]"
                ")=C([H])C([H])=C3C([H])=N->%13C3=C([H])C([H])=C(C([H]"
                ")=C3[H])C([H])(C3=C([H])C([H])=C(C([H])=C3[H])N->9=C("
                "[H])C3=C([H])C([H])=C([H])C([H])=N->%103)C([H])(C3=C("
                "[H])C([H])=C(C([H])=C3[H])N->4=C([H])C3=C([H])C([H])="
                "C([H])C([H])=N->83)C3=C([H])C([H])=C(C([H])=C3[H])N->"
                "6=C([H])C3=C([H])C([H])=C([H])C([H])=N->%113)<-N3=C(["
                "H])C([H])=C([H])C([H])=C3C([H])=N->%14C3=C([H])C([H])"
                "=C1C([H])=C3[H])<-N1=C([H])C([H])=C([H])C([H])=C1C=2["
                "H])<-N1=C([H])C([H])=C([H])C([H])=C1C=7[H]"
            ),
            name=name,
        ),
    ),
)
def metal_cage_m6l2l3_prism(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
