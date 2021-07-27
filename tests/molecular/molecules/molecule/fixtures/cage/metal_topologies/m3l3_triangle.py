import pytest
import stk

from .building_blocks import get_metal_atom
from ....case_data import CaseData


def _get_palladium_bi_1() -> stk.BuildingBlock:
    return stk.BuildingBlock(
        smiles='[H]N([H])C([H])([H])C([H])([H])N([H])[H]',
        functional_groups=[
            stk.SmartsFunctionalGroupFactory(
                smarts='[#7]~[#6]',
                bonders=(0, ),
                deleters=(),
            ),
        ],
    )


def _get_palladium_cispbi_sqpl() -> stk.BuildingBlock:
    molecule = stk.ConstructedMolecule(
        topology_graph=stk.metal_complex.CisProtectedSquarePlanar(
            metals={get_metal_atom(): 0},
            ligands={_get_palladium_bi_1(): 0},
            reaction_factory=stk.DativeReactionFactory(
                reaction_factory=stk.GenericReactionFactory(
                    bond_orders={
                        frozenset({
                            stk.GenericFunctionalGroup,
                            stk.SingleAtom
                        }): 9,
                    },
                ),
            ),
        ),
    )

    return stk.BuildingBlock.init_from_molecule(
        molecule=molecule,
        functional_groups=[
            stk.SmartsFunctionalGroupFactory(
                smarts='[Pd]~[#7]',
                bonders=(0, ),
                deleters=(),
                placers=(0, 1),
            ),
        ],
    )


def _get_linker() -> stk.BuildingBlock:
    return stk.BuildingBlock(
        smiles=(
            '[H]C1=NC([H])=C([H])C(C2=C([H])C([H])=NC([H])=C2[H])'
            '=C1[H]'
        ),
        functional_groups=[
            stk.SmartsFunctionalGroupFactory(
                smarts='[#6]~[#7X2]~[#6]',
                bonders=(1, ),
                deleters=(),
            ),
        ],
    )


@pytest.fixture(
    scope='session',
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.M3L3Triangle(
                    corners=_get_palladium_cispbi_sqpl(),
                    linkers=_get_linker(),
                    reaction_factory=stk.DativeReactionFactory(
                        reaction_factory=stk.GenericReactionFactory(
                            bond_orders={
                                frozenset({
                                    stk.GenericFunctionalGroup,
                                    stk.GenericFunctionalGroup
                                }): 9,
                            },
                        ),
                    ),
                ),
            ),
            smiles=(
                '[H]C1=C2C([H])=C([H])N(->[Pd+2]3(<-N4=C([H])C([H])=C'
                '(C([H])=C4[H])C4=C([H])C([H])=N(->[Pd+2]5(<-N6=C([H]'
                ')C([H])=C(C([H])=C6[H])C6=C([H])C([H])=N(->[Pd+2]7(<'
                '-N8=C([H])C([H])=C2C([H])=C8[H])<-N([H])([H])C([H])('
                '[H])C([H])([H])N->7([H])[H])C([H])=C6[H])<-N([H])([H'
                '])C([H])([H])C([H])([H])N->5([H])[H])C([H])=C4[H])<-'
                'N([H])([H])C([H])([H])C([H])([H])N->3([H])[H])=C1[H]'
            ),
            name=name,
        ),
    ),
)
def metal_cage_m3l3_triangle(request) -> CaseData:
    return request.param(request.node.originalname)
