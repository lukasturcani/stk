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
                "[H]C1=C([H])C2=C([H])C(=C1[H])C1=C([H])C([H])=N(->["
                "Pd+2]34<-N5=C([H])C([H])=C(C([H])=C5[H])C5=C([H])C("
                "[H])=C([H])C(=C5[H])C5=C([H])C([H])=N(->[Pd+2]67<-N"
                "8=C([H])C([H])=C(C([H])=C8[H])C8=C([H])C([H])=C([H]"
                ")C(=C8[H])C8=C([H])C([H])=N(->[Pd+2]9(<-N%10=C([H])"
                "C([H])=C(C([H])=C%10[H])C%10=C([H])C(=C([H])C([H])="
                "C%10[H])C%10=C([H])C([H])=N->3C([H])=C%10[H])<-N3=C"
                "([H])C([H])=C(C([H])=C3[H])C3=C([H])C(=C([H])C([H])"
                "=C3[H])C3=C([H])C([H])=N(->[Pd+2]%10(<-N%11=C([H])C"
                "([H])=C(C([H])=C%11[H])C%11=C([H])C([H])=C([H])C(=C"
                "%11[H])C%11=C([H])C([H])=N(->[Pd+2](<-N%12=C([H])C("
                "[H])=C(C([H])=C%12[H])C%12=C([H])C([H])=C([H])C(=C%"
                "12[H])C%12=C([H])C([H])=N->6C([H])=C%12[H])(<-N6=C("
                "[H])C([H])=C(C([H])=C6[H])C6=C([H])C([H])=C([H])C(="
                "C6[H])C6=C([H])C([H])=N->9C([H])=C6[H])<-N6=C([H])C"
                "([H])=C(C([H])=C6[H])C6=C([H])C([H])=C([H])C(=C6[H]"
                ")C6=C([H])C([H])=N(->[Pd+2](<-N9=C([H])C([H])=C(C(["
                "H])=C9[H])C9=C([H])C(=C([H])C([H])=C9[H])C9=C([H])C"
                "([H])=N->4C([H])=C9[H])(<-N4=C([H])C([H])=C(C([H])="
                "C4[H])C4=C([H])C(=C([H])C([H])=C4[H])C4=C([H])C([H]"
                ")=N->%10C([H])=C4[H])<-N4=C([H])C([H])=C(C([H])=C4["
                "H])C4=C([H])C(=C([H])C([H])=C4[H])C4=C([H])C([H])=N"
                "->7C([H])=C4[H])C([H])=C6[H])C([H])=C%11[H])<-N4=C("
                "[H])C([H])=C2C([H])=C4[H])C([H])=C3[H])C([H])=C8[H]"
                ")C([H])=C5[H])C([H])=C1[H]"
            ),
            name=name,
        ),
    ),
)
def metal_cage_m6l12_cube(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
