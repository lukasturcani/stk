import pytest
import numpy as np
import stk

from .case_data import CaseData

bb1 = stk.BuildingBlock('[C+2][N+]Br', [stk.BromoFactory()])
bb2 = stk.BuildingBlock('IS[O+]', [stk.IodoFactory()])


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.BuildingBlock(
                smiles='Br[C+2][N+]Cl',
                functional_groups=[stk.BromoFactory()],
                placer_ids=(0, ),
            ),
            result=stk.BuildingBlock.init(
                atoms=(
                    stk.Cl(0),
                    stk.Br(1),
                    stk.C(2, 2),
                    stk.N(3, 1),
                ),
                bonds=(
                    stk.Bond(stk.Br(1), stk.C(2, 2), 1),
                    stk.Bond(stk.C(2, 2), stk.N(3, 1), 1),
                    stk.Bond(stk.N(3, 1), stk.Cl(0), 1),
                ),
                position_matrix=np.array([
                ]),
                functional_groups=(
                    stk.Bromo(
                        bromine=stk.Br(1),
                        atom=stk.C(2, 2),
                        bonders=(stk.C(2, 2), ),
                        deleters=(stk.Br(1), ),
                    ),
                ),
                placer_ids=(1, ),
            )
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(bb1, bb2),
                    repeating_unit='AB',
                    num_repeating_units=1,
                ),
            ),
            result=stk.ConstructedMolecule.init(
                atoms=(
                    stk.C(0, 2),
                    stk.O(1, 1),
                    stk.N(2, 1),
                    stk.S(3),
                ),
                bonds=(
                    stk.Bond(stk.C(0, 2), stk.N(2, 1), 1),
                    stk.Bond(stk.S(3), stk.O(1, 1), 1),
                    stk.Bond(stk.N(2, 1), stk.S(3), 1),
                ),
                position_matrix=np.array([]),
                atom_infos=(
                    stk.AtomInfo(
                        atom=stk.C(0, 2),
                        building_block=bb1,
                        building_block_id=0,
                    ),
                    stk.AtomInfo(
                        atom=stk.O(1, 1),
                        building_block=bb2,
                        building_block_id=1,
                    ),
                    stk.AtomInfo(
                        atom=stk.N(2, 1),
                        building_block=bb1,
                        building_block_id=0,
                    ),
                    stk.AtomInfo(
                        atom=stk.S(3),
                        building_block=bb2,
                        building_block_id=1,
                    ),
                ),
                bond_infos=(
                    stk.BondInfo(
                        bond=stk.Bond(stk.C(0, 2), stk.N(2, 1), 1),
                        building_block=bb1,
                        building_block_id=0,
                    ),
                    stk.BondInfo(
                        bond=stk.Bond(stk.S(3), stk.O(1, 1), 1),
                        building_block=bb2,
                        building_block_id=1,
                    ),
                    stk.BondInfo(
                        bond=stk.Bond(stk.N(2, 1), stk.S(3), 1),
                        building_block=None,
                        building_block_id=None,
                    ),
                ),
                num_building_blocks={
                    bb1: 1,
                    bb2: 1,
                },
            ),
        ),
    ),
)
def case_data(request):
    return request.param
