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
                stk.cage.M4L4Tetrahedron(
                    building_blocks={
                        iron_complex: range(4),
                        tritopic_linker: range(4, 8)
                    },
                )
            ),
            smiles=(
                '[H]C1=C([H])C([H])=N2->[Fe+2]3456<-N7=C([H])C([H])=C'
                '([H])C([H])=C7C([H])=N->3C3=C([H])C([H])=C(C([H])=C3'
                '[H])N3C7=C([H])C([H])=C(C([H])=C7[H])N7->[Fe+2]89%10'
                '(<-N%11=C([H])C([H])=C([H])C([H])=C%11C([H])=N->8C8='
                'C([H])C([H])=C(C([H])=C8[H])N(C8=C([H])C([H])=C(C([H'
                '])=C8[H])N->4=C([H])C2=C1[H])C1=C([H])C([H])=C(C([H]'
                ')=C1[H])N1->[Fe+2]248(<-N%11=C([H])C([H])=C([H])C([H'
                '])=C%11C=1[H])<-N1=C([H])C([H])=C([H])C([H])=C1C([H]'
                ')=N->2C1=C([H])C([H])=C(C([H])=C1[H])N(C1=C([H])C([H'
                '])=C(C([H])=C1[H])N->5=C([H])C1=C([H])C([H])=C([H])C'
                '([H])=N->61)C1=C([H])C([H])=C(C([H])=C1[H])N1->[Fe+2'
                ']25(<-N6=C([H])C([H])=C([H])C([H])=C6C([H])=N->2C2=C'
                '([H])C([H])=C3C([H])=C2[H])(<-N2=C([H])C([H])=C([H])'
                'C([H])=C2C=1[H])<-N1=C([H])C([H])=C([H])C([H])=C1C(['
                'H])=N->5C1=C([H])C([H])=C(C([H])=C1[H])N(C1=C([H])C('
                '[H])=C(C([H])=C1[H])N->9=C([H])C1=C([H])C([H])=C([H]'
                ')C([H])=N->%101)C1=C([H])C([H])=C(C([H])=C1[H])N->4='
                'C([H])C1=C([H])C([H])=C([H])C([H])=N->81)<-N1=C([H])'
                'C([H])=C([H])C([H])=C1C=7[H]'
            ),
        ),
    ),
)
def mcage_m4l4_tetrahedron(request):
    return request.param
