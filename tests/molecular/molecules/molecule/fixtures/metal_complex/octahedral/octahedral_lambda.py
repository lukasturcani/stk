import pytest
import stk

from ....case_data import CaseData
from ...building_blocks import get_fe_atom, get_iron_bi_1


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                stk.metal_complex.OctahedralLambda(
                    metals={get_fe_atom(): 0},
                    ligands={get_iron_bi_1(): (0, 1, 2)},
                    reaction_factory=stk.DativeReactionFactory(
                        stk.GenericReactionFactory(
                            bond_orders={
                                frozenset(
                                    {
                                        stk.GenericFunctionalGroup,
                                        stk.SingleAtom,
                                    }
                                ): 9
                            }
                        )
                    ),
                )
            ),
            smiles=(
                "[H]C1=C([H])C([H])=N2->[Fe+2]34(<-N(Br)=C([H])C2=C1[H"
                "])(<-N(Br)=C([H])C1=C([H])C([H])=C([H])C([H])=N->31)"
                "<-N(Br)=C([H])C1=C([H])C([H])=C([H])C([H])=N->41"
            ),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                stk.metal_complex.OctahedralLambda(
                    metals=get_fe_atom(),
                    ligands=get_iron_bi_1(),
                    reaction_factory=stk.DativeReactionFactory(
                        stk.GenericReactionFactory(
                            bond_orders={
                                frozenset(
                                    {
                                        stk.GenericFunctionalGroup,
                                        stk.SingleAtom,
                                    }
                                ): 9
                            }
                        )
                    ),
                )
            ),
            smiles=(
                "[H]C1=C([H])C([H])=N2->[Fe+2]34(<-N(Br)=C([H])C2=C1[H"
                "])(<-N(Br)=C([H])C1=C([H])C([H])=C([H])C([H])=N->31)"
                "<-N(Br)=C([H])C1=C([H])C([H])=C([H])C([H])=N->41"
            ),
            name=name,
        ),
    ),
)
def metal_complex_octahedral_lambda(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
