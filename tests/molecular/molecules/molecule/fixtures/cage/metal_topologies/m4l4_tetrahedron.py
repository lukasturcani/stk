import pytest
import stk

from ....case_data import CaseData
from ...building_blocks import get_iron_complex, get_tritopic_linker


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.M4L4Tetrahedron(
                    building_blocks={
                        get_iron_complex(): range(4),
                        get_tritopic_linker(): range(4, 8),
                    },
                ),
            ),
            smiles=(
                "[H]C1=C([H])C([H])=N2->[Fe+2]3456<-N7=C([H])C([H])=C"
                "([H])C([H])=C7C([H])=N->3C3=C([H])C([H])=C(C([H])=C3"
                "[H])N3C7=C([H])C([H])=C(C([H])=C7[H])N7->[Fe+2]89%10"
                "(<-N%11=C([H])C([H])=C([H])C([H])=C%11C([H])=N->8C8="
                "C([H])C([H])=C(C([H])=C8[H])N(C8=C([H])C([H])=C(C([H"
                "])=C8[H])N->4=C([H])C2=C1[H])C1=C([H])C([H])=C(C([H]"
                ")=C1[H])N1->[Fe+2]248(<-N%11=C([H])C([H])=C([H])C([H"
                "])=C%11C=1[H])<-N1=C([H])C([H])=C([H])C([H])=C1C([H]"
                ")=N->2C1=C([H])C([H])=C(C([H])=C1[H])N(C1=C([H])C([H"
                "])=C(C([H])=C1[H])N->5=C([H])C1=C([H])C([H])=C([H])C"
                "([H])=N->61)C1=C([H])C([H])=C(C([H])=C1[H])N1->[Fe+2"
                "]25(<-N6=C([H])C([H])=C([H])C([H])=C6C([H])=N->2C2=C"
                "([H])C([H])=C3C([H])=C2[H])(<-N2=C([H])C([H])=C([H])"
                "C([H])=C2C=1[H])<-N1=C([H])C([H])=C([H])C([H])=C1C(["
                "H])=N->5C1=C([H])C([H])=C(C([H])=C1[H])N(C1=C([H])C("
                "[H])=C(C([H])=C1[H])N->9=C([H])C1=C([H])C([H])=C([H]"
                ")C([H])=N->%101)C1=C([H])C([H])=C(C([H])=C1[H])N->4="
                "C([H])C1=C([H])C([H])=C([H])C([H])=N->81)<-N1=C([H])"
                "C([H])=C([H])C([H])=C1C=7[H]"
            ),
            name=name,
        ),
    ),
)
def metal_cage_m4l4_tetrahedron(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
