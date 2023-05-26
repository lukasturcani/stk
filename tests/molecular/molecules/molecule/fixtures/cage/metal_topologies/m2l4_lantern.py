import pytest
import stk

from ....case_data import CaseData
from ...building_blocks import get_linker, get_pd_atom


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                stk.cage.M2L4Lantern(
                    building_blocks={
                        get_pd_atom(): range(2),
                        get_linker(): range(2, 6),
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
                "[H]C1=C([H])C2=C([H])C(=C1[H])C1=C([H])C([H])=N(->["
                "Pd+2]34<-N5=C([H])C([H])=C(C([H])=C5[H])C5=C([H])C("
                "[H])=C([H])C(=C5[H])C5=C([H])C([H])=N(->[Pd+2](<-N6"
                "=C([H])C([H])=C2C([H])=C6[H])(<-N2=C([H])C([H])=C(C"
                "([H])=C2[H])C2=C([H])C(=C([H])C([H])=C2[H])C2=C([H]"
                ")C([H])=N->3C([H])=C2[H])<-N2=C([H])C([H])=C(C([H])"
                "=C2[H])C2=C([H])C(=C([H])C([H])=C2[H])C2=C([H])C([H"
                "])=N->4C([H])=C2[H])C([H])=C5[H])C([H])=C1[H]"
            ),
            name=name,
        ),
    ),
)
def metal_cage_m2l4_lantern(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
