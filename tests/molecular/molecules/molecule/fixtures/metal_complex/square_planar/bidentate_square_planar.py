import pytest
import stk
from rdkit.Chem import AllChem as rdkit

from ....case_data import CaseData


atom = rdkit.MolFromSmiles('[Pd+2]')
atom.AddConformer(rdkit.Conformer(atom.GetNumAtoms()))

_palladium_atom = stk.BuildingBlock.init_from_rdkit_mol(atom)
atom_0, = _palladium_atom.get_atoms(0)
_palladium_atom = _palladium_atom.with_functional_groups(
    (stk.SingleAtom(atom_0) for i in range(4))
)

_bi_1 = stk.BuildingBlock(
    smiles='NCCN',
    functional_groups=[
        stk.SmartsFunctionalGroupFactory(
            smarts='[#7]~[#6]',
            bonders=(0, ),
            deleters=(),
        ),
    ]
)


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.ConstructedMolecule(
                stk.metal_complex.BidentateSquarePlanar(
                    metals={_palladium_atom: 0},
                    ligands={_bi_1: (0, 1)},
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
                '[H]C1([H])C([H])([H])N([H])([H])->[Pd+2]2(<-N1([H])'
                '[H])<-N([H])([H])C([H])([H])C([H])([H])N->2([H])[H]'
            ),
        ),

        CaseData(
            molecule=stk.ConstructedMolecule(
                stk.metal_complex.BidentateSquarePlanar(
                    metals=_palladium_atom,
                    ligands=_bi_1,
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
                '[H]C1([H])C([H])([H])N([H])([H])->[Pd+2]2(<-N1([H])'
                '[H])<-N([H])([H])C([H])([H])C([H])([H])N->2([H])[H]'
            ),
        ),
    ),
)
def metal_complex_bidentatesquareplanar(request):
    return request.param
