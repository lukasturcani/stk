import pytest
import stk

from ....case_data import CaseData


metal_atom = stk.BuildingBlock(
    smiles='[Fe+2]',
    functional_groups=(
        stk.SingleAtom(stk.Fe(0, charge=2))
        for i in range(6)
    ),
    position_matrix=([0, 0, 0], ),
)

tritopic_linker = stk.BuildingBlock(
    smiles=(
        '[H]C1=C([H])C(N(C2=C([H])C([H])=C(Br)C([H])=C2[H])C2=C([H])C('
        '[H])=C(Br)C([H])=C2[H])=C([H])C([H])=C1Br'
    ),
    functional_groups=[stk.BromoFactory()]
)
tetratopic_linker = stk.BuildingBlock(
    smiles=(
        '[H]C1=C([H])C(C(C2=C([H])C([H])=C(Br)C([H])=C2[H])C(C2=C([H])'
        'C([H])=C(Br)C([H])=C2[H])C2=C([H])C([H])=C(Br)C([H])=C2[H])=C'
        '([H])C([H])=C1Br'
    ),
    functional_groups=[stk.BromoFactory()]
)
complex_ligand = stk.BuildingBlock(
    smiles='[H]C1=NC(C([H])=NBr)=C([H])C([H])=C1[H]',
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
iron_complex = stk.ConstructedMolecule(
    stk.metal_complex.OctahedralDelta(
        metals={metal_atom: 0},
        ligands={complex_ligand: (0, 1, 2)},
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
)
iron_complex = stk.BuildingBlock.init_from_molecule(
    molecule=iron_complex,
    functional_groups=[stk.BromoFactory()]
)


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.ConstructedMolecule(
                stk.cage.M6L2L3Prism(
                    building_blocks={
                        iron_complex: range(6),
                        tritopic_linker: range(6, 8),
                        tetratopic_linker: range(8, 11),
                    },
                )
            ),
            smiles=(
                '[H]C1=C([H])C([H])=N2->[Fe+2]3456<-N7=C([H])C([H])=C('
                '[H])C([H])=C7C([H])=N->3C3=C([H])C([H])=C(C([H])=C3[H'
                '])C3([H])C7=C([H])C([H])=C(C([H])=C7[H])N7->[Fe+2]89%'
                '10(<-N%11=C([H])C([H])=C([H])C([H])=C%11C([H])=N->8C8'
                '=C([H])C([H])=C(C([H])=C8[H])N(C8=C([H])C([H])=C(C([H'
                '])=C8[H])N->4=C([H])C2=C1[H])C1=C([H])C([H])=C(C([H])'
                '=C1[H])N1->[Fe+2]248(<-N%11=C([H])C([H])=C([H])C([H])'
                '=C%11C=1[H])<-N1=C([H])C([H])=C([H])C([H])=C1C([H])=N'
                '->2C1=C([H])C([H])=C(C([H])=C1[H])C([H])(C1=C([H])C(['
                'H])=C(C([H])=C1[H])N->5=C([H])C1=C([H])C([H])=C([H])C'
                '([H])=N->61)C1([H])C2=C([H])C([H])=C(C([H])=C2[H])N2-'
                '>[Fe+2]56%11(<-N%12=C([H])C([H])=C([H])C([H])=C%12C(['
                'H])=N->5C5=C([H])C([H])=C(C([H])=C5[H])N5C%12=C([H])C'
                '([H])=C(C([H])=C%12[H])N%12->[Fe+2]%13%14(<-N%15=C([H'
                '])C([H])=C([H])C([H])=C%15C=%12[H])(<-N%12=C([H])C([H'
                '])=C([H])C([H])=C%12C([H])=N->%13C%12=C([H])C([H])=C('
                'C([H])=C%12[H])C3([H])C3=C([H])C([H])=C(C([H])=C3[H])'
                'N3->[Fe+2]%12%13(<-N%15=C([H])C([H])=C([H])C([H])=C%1'
                '5C([H])=N->%12C%12=C([H])C([H])=C5C([H])=C%12[H])(<-N'
                '5=C([H])C([H])=C([H])C([H])=C5C=3[H])<-N3=C([H])C([H]'
                ')=C([H])C([H])=C3C([H])=N->%13C3=C([H])C([H])=C(C([H]'
                ')=C3[H])C([H])(C3=C([H])C([H])=C(C([H])=C3[H])N->9=C('
                '[H])C3=C([H])C([H])=C([H])C([H])=N->%103)C([H])(C3=C('
                '[H])C([H])=C(C([H])=C3[H])N->4=C([H])C3=C([H])C([H])='
                'C([H])C([H])=N->83)C3=C([H])C([H])=C(C([H])=C3[H])N->'
                '6=C([H])C3=C([H])C([H])=C([H])C([H])=N->%113)<-N3=C(['
                'H])C([H])=C([H])C([H])=C3C([H])=N->%14C3=C([H])C([H])'
                '=C1C([H])=C3[H])<-N1=C([H])C([H])=C([H])C([H])=C1C=2['
                'H])<-N1=C([H])C([H])=C([H])C([H])=C1C=7[H]'
            ),
        ),
    ),
)
def mcage_m6l2l3_prism(request):
    return request.param
