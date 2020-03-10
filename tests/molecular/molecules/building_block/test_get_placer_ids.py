import pytest
import rdkit.Chem.AllChem as rdkit
import numpy as np
import stk
import itertools as it


class CaseData:
    def __init__(self, building_block, placer_ids):
        self.building_block = building_block
        self.placer_ids = placer_ids


@pytest.fixture(
    params=(
        CaseData(
            building_block=stk.BuildingBlock('Br[C+2][C+2]Br'),
            placer_ids=(0, 1, 2, 3),
        ),
        CaseData(
            building_block=stk.BuildingBlock(
                smiles='BrCCCBr',
                functional_groups=[stk.BromoFactory()],
                placer_ids=(),
            ),
            placer_ids=(),
        ),
        CaseData(
            building_block=stk.BuildingBlock(
                smiles='BrCCCBr',
                functional_groups=[stk.BromoFactory()],
                placer_ids=(2, ),
            ),
            placer_ids=(2, ),
        ),
        CaseData(
            building_block=stk.BuildingBlock(
                smiles='Br[C+2][C+2]Br',
                functional_groups=[stk.IodoFactory()],
            ),
            placer_ids=(0, 1, 2, 3),
        ),
        CaseData(
            building_block=stk.BuildingBlock(
                smiles='Br[C+2][C+2]Br',
                functional_groups=[stk.BromoFactory()],
            ),
            placer_ids=(1, 2),
        ),
        CaseData(
            building_block=stk.BuildingBlock.init(
                atoms=(stk.Br(0), stk.C(1), stk.C(2), stk.Br(3)),
                bonds=(
                    stk.Bond(stk.Br(0), stk.C(1), 1),
                    stk.Bond(stk.C(1), stk.C(2), 1),
                    stk.Bond(stk.C(2), stk.Br(3), 1),
                ),
                position_matrix=np.array([
                    [0., 0., 0.],
                    [2., 0., 0.],
                    [3., 0., 0.],
                    [4., 0., 0.],
                ]),
                functional_groups=(
                    stk.Bromo(stk.Br(0), stk.C(1), (stk.C(1), ), ()),
                    stk.Bromo(stk.Br(3), stk.C(2), (stk.C(2), ), ()),
                ),
                placer_ids=(),
            ),
            placer_ids=(),
        ),
        CaseData(
            building_block=stk.BuildingBlock.init(
                atoms=(stk.Br(0), stk.C(1), stk.C(2), stk.Br(3)),
                bonds=(
                    stk.Bond(stk.Br(0), stk.C(1), 1),
                    stk.Bond(stk.C(1), stk.C(2), 1),
                    stk.Bond(stk.C(2), stk.Br(3), 1),
                ),
                position_matrix=np.array([
                    [0., 0., 0.],
                    [2., 0., 0.],
                    [3., 0., 0.],
                    [4., 0., 0.],
                ]),
                functional_groups=(
                    stk.Bromo(stk.Br(0), stk.C(1), (stk.C(1), ), ()),
                    stk.Bromo(stk.Br(3), stk.C(2), (stk.C(2), ), ()),
                ),
                placer_ids=(3, ),
            ),
            placer_ids=(3, ),
        ),
        CaseData(
            building_block=stk.BuildingBlock.init(
                atoms=(stk.Br(0), stk.C(1), stk.C(2), stk.Br(3)),
                bonds=(
                    stk.Bond(stk.Br(0), stk.C(1), 1),
                    stk.Bond(stk.C(1), stk.C(2), 1),
                    stk.Bond(stk.C(2), stk.Br(3), 1),
                ),
                position_matrix=np.array([
                    [0., 0., 0.],
                    [2., 0., 0.],
                    [3., 0., 0.],
                    [4., 0., 0.],
                ]),
                functional_groups=[stk.IodoFactory()],
            ),
            placer_ids=(0, 1, 2, 3),
        ),
        CaseData(
            building_block=stk.BuildingBlock.init(
                atoms=(stk.Br(0), stk.C(1), stk.C(2), stk.Br(3)),
                bonds=(
                    stk.Bond(stk.Br(0), stk.C(1), 1),
                    stk.Bond(stk.C(1), stk.C(2), 1),
                    stk.Bond(stk.C(2), stk.Br(3), 1),
                ),
                position_matrix=np.array([
                    [0., 0., 0.],
                    [2., 0., 0.],
                    [3., 0., 0.],
                    [4., 0., 0.],
                ]),
                functional_groups=(
                    stk.Bromo(stk.Br(0), stk.C(1), (stk.C(1), ), ()),
                    stk.Bromo(stk.Br(3), stk.C(2), (stk.C(2), ), ()),
                ),
            ),
            placer_ids=(1, 2),
        ),
        CaseData(
            building_block=stk.BuildingBlock.init_from_molecule(
                molecule=stk.BuildingBlock('Br[C+2][C+2]Br'),
                functional_groups=[stk.IodoFactory()],
            ),
            placer_ids=(0, 1, 2, 3),
        ),
        CaseData(
            building_block=stk.BuildingBlock.init_from_molecule(
                molecule=stk.BuildingBlock('Br[C+2][C+2]Br'),
                functional_groups=[stk.BromoFactory()],
            ),
            placer_ids=(1, 2),
        ),
        CaseData(
            building_block=stk.BuildingBlock.init_from_molecule(
                molecule=stk.BuildingBlock('Br[C+2][C+2]Br'),
                functional_groups=[stk.BromoFactory()],
                placer_ids=(3, ),
            ),
            placer_ids=(3, ),
        ),
        CaseData(
            building_block=stk.BuildingBlock.init_from_rdkit_mol(
                molecule=(
                    stk.BuildingBlock('Br[C+2][C+2]Br').to_rdkit_mol()
                ),
                functional_groups=[stk.IodoFactory()],
            ),
            placer_ids=(0, 1, 2, 3),
        ),
        CaseData(
            building_block=stk.BuildingBlock.init_from_rdkit_mol(
                molecule=(
                    stk.BuildingBlock('Br[C+2][C+2]Br').to_rdkit_mol()
                ),
                functional_groups=[stk.BromoFactory()],
            ),
            placer_ids=(1, 2),
        ),
        CaseData(
            building_block=stk.BuildingBlock.init_from_rdkit_mol(
                molecule=(
                    stk.BuildingBlock('Br[C+2][C+2]Br').to_rdkit_mol()
                ),
                functional_groups=[stk.BromoFactory()],
                placer_ids=(3, ),
            ),
            placer_ids=(3, ),
        ),
    ),
)
def case_data(request):
    return request.param


def test_get_placer_ids(case_data):
    _test_get_placer_ids(
        building_block=case_data.building_block,
        placer_ids=case_data.placer_ids,
    )


def _test_get_placer_ids(building_block, placer_ids):
    for placer1, placer2 in it.zip_longest(
        building_block.get_placer_ids(),
        placer_ids,
    ):
        assert placer1 == placer2
