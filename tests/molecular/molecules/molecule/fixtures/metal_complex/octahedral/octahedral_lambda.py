import pytest
import stk
from rdkit.Chem import AllChem as rdkit

from ....case_data import CaseData


atom = rdkit.MolFromSmiles('[Fe+2]')
atom.AddConformer(rdkit.Conformer(atom.GetNumAtoms()))

_iron_atom = stk.BuildingBlock.init_from_rdkit_mol(atom)
atom_0, = _iron_atom.get_atoms(0)
_iron_atom = _iron_atom.with_functional_groups(
    (stk.SingleAtom(atom_0) for i in range(6))
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
                stk.metal_complex.OctahedralLambda(
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
                stk.metal_complex.OctahedralLambda(
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
def metal_complex_octahedrallambda(request):
    return request.param
