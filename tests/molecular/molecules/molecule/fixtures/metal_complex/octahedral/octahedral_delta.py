import pytest
import stk

from ....case_data import CaseData


_iron_atom = stk.BuildingBlock(
    smiles='[Fe+2]',
    functional_groups=(
        stk.SingleAtom(stk.Fe(0, charge=2))
        for i in range(6)
    ),
    position_matrix=([0, 0, 0], ),
)

_iron_bi_1 = stk.BuildingBlock(
    smiles='BrN=Cc1ccccn1',
    functional_groups=[
        stk.SmartsFunctionalGroupFactory(
            smarts='[#6]~[#7X2]~[#35]',
            bonders=(1, ),
            deleters=(),
        ),
        stk.SmartsFunctionalGroupFactory(
            smarts='[#6]~[#7X2]~[#6]',
            bonders=(1, ),
            deleters=(),
        ),
    ]
)


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.ConstructedMolecule(
                stk.metal_complex.OctahedralDelta(
                    metals={_iron_atom: 0},
                    ligands={_iron_bi_1: (0, 1, 2)},
                    reaction_factory=stk.DativeReactionFactory(
                        stk.GenericReactionFactory(
                            bond_orders={
                                frozenset({
                                    stk.GenericFunctionalGroup,
                                    stk.SingleAtom
                                }): 9
                            }
                        )
                    )
                )
            ),
            smiles=(
                '[H]C1=C([H])C([H])=N2->[Fe+2]34(<-N(Br)=C([H])C2=C1[H'
                '])(<-N(Br)=C([H])C1=C([H])C([H])=C([H])C([H])=N->31)'
                '<-N(Br)=C([H])C1=C([H])C([H])=C([H])C([H])=N->41'
            ),
        ),

        CaseData(
            molecule=stk.ConstructedMolecule(
                stk.metal_complex.OctahedralDelta(
                    metals=_iron_atom,
                    ligands=_iron_bi_1,
                    reaction_factory=stk.DativeReactionFactory(
                        stk.GenericReactionFactory(
                            bond_orders={
                                frozenset({
                                    stk.GenericFunctionalGroup,
                                    stk.SingleAtom
                                }): 9
                            }
                        )
                    )
                )
            ),
            smiles=(
                '[H]C1=C([H])C([H])=N2->[Fe+2]34(<-N(Br)=C([H])C2=C1[H'
                '])(<-N(Br)=C([H])C1=C([H])C([H])=C([H])C([H])=N->31)'
                '<-N(Br)=C([H])C1=C([H])C([H])=C([H])C([H])=N->41'
            ),
        ),
    ),
)
def metal_complex_octahedraldelta(request):
    return request.param
