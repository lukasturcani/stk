import pytest
import stk

from ....case_data import CaseData
from ...building_blocks import get_fe_atom, get_mo_1


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.metal_complex.Octahedral(
                    metals={get_fe_atom(): 0},
                    ligands={get_mo_1(): (0, 1, 2, 3, 4, 5)},
                    reaction_factory=stk.DativeReactionFactory(
                        stk.GenericReactionFactory(
                            bond_orders={
                                frozenset(
                                    {
                                        stk.GenericFunctionalGroup,
                                        stk.SingleAtom,
                                    }
                                ): 9,
                            }
                        )
                    ),
                )
            ),
            smiles=(
                "[H]C1=C([H])N(->[Fe+2](<-N2=C([H])C3=C(C([H])=C2[H])C"
                "([H])([H])C([H])([H])C([H])([H])C([H])([H])C3([H])[H]"
                ")(<-N2=C([H])C3=C(C([H])=C2[H])C([H])([H])C([H])([H])"
                "C([H])([H])C([H])([H])C3([H])[H])(<-N2=C([H])C3=C(C("
                "[H])=C2[H])C([H])([H])C([H])([H])C([H])([H])C([H])(["
                "H])C3([H])[H])(<-N2=C([H])C3=C(C([H])=C2[H])C([H])(["
                "H])C([H])([H])C([H])([H])C([H])([H])C3([H])[H])<-N2="
                "C([H])C3=C(C([H])=C2[H])C([H])([H])C([H])([H])C([H])"
                "([H])C([H])([H])C3([H])[H])=C([H])C2=C1C([H])([H])C("
                "[H])([H])C([H])([H])C([H])([H])C2([H])[H]"
            ),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.metal_complex.Octahedral(
                    metals=get_fe_atom(),
                    ligands=get_mo_1(),
                    reaction_factory=stk.DativeReactionFactory(
                        reaction_factory=stk.GenericReactionFactory(
                            bond_orders={
                                frozenset(
                                    {
                                        stk.GenericFunctionalGroup,
                                        stk.SingleAtom,
                                    }
                                ): 9,
                            }
                        )
                    ),
                )
            ),
            smiles=(
                "[H]C1=C([H])N(->[Fe+2](<-N2=C([H])C3=C(C([H])=C2[H])C"
                "([H])([H])C([H])([H])C([H])([H])C([H])([H])C3([H])[H]"
                ")(<-N2=C([H])C3=C(C([H])=C2[H])C([H])([H])C([H])([H])"
                "C([H])([H])C([H])([H])C3([H])[H])(<-N2=C([H])C3=C(C("
                "[H])=C2[H])C([H])([H])C([H])([H])C([H])([H])C([H])(["
                "H])C3([H])[H])(<-N2=C([H])C3=C(C([H])=C2[H])C([H])(["
                "H])C([H])([H])C([H])([H])C([H])([H])C3([H])[H])<-N2="
                "C([H])C3=C(C([H])=C2[H])C([H])([H])C([H])([H])C([H])"
                "([H])C([H])([H])C3([H])[H])=C([H])C2=C1C([H])([H])C("
                "[H])([H])C([H])([H])C([H])([H])C2([H])[H]"
            ),
            name=name,
        ),
    ),
)
def metal_complex_octahedral(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
