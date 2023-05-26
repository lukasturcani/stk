import pytest
import stk

from ....case_data import CaseData
from ...building_blocks import get_linker, get_pd_atom


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.M3L6(
                    building_blocks={
                        get_pd_atom(): range(3),
                        get_linker(): range(3, 9),
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
                "[H])=C([H])C(=C5[H])C5=C([H])C([H])=N(->[Pd+2]6(<-N"
                "7=C([H])C([H])=C(C([H])=C7[H])C7=C([H])C([H])=C([H]"
                ")C(=C7[H])C7=C([H])C([H])=N(->[Pd+2](<-N8=C([H])C(["
                "H])=C(C([H])=C8[H])C8=C([H])C(=C([H])C([H])=C8[H])C"
                "8=C([H])C([H])=N->3C([H])=C8[H])(<-N3=C([H])C([H])="
                "C(C([H])=C3[H])C3=C([H])C(=C([H])C([H])=C3[H])C3=C("
                "[H])C([H])=N->4C([H])=C3[H])<-N3=C([H])C([H])=C(C(["
                "H])=C3[H])C3=C([H])C(=C([H])C([H])=C3[H])C3=C([H])C"
                "([H])=N->6C([H])=C3[H])C([H])=C7[H])<-N3=C([H])C([H"
                "])=C2C([H])=C3[H])C([H])=C5[H])C([H])=C1[H]"
            ),
            name=name,
        ),
    ),
)
def metal_cage_m3l6(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
