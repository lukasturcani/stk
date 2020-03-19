import pytest
import stk

from ..case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            building_block=stk.BuildingBlock('Br[C+2][C+2]Br'),
            functional_groups=(),
            core_atoms=(stk.Br(0), stk.C(1), stk.C(2), stk.Br(3)),
            placers=(stk.Br(0), stk.C(1), stk.C(2), stk.Br(3)),
        ),
        CaseData(
            building_block=stk.BuildingBlock(
                smiles='Br[C+2][C+2]Br',
                functional_groups=[stk.BromoFactory()],
            ),
            functional_groups=(
                stk.Bromo(
                    bromine=stk.Br(0),
                    atom=stk.C(1),
                    bonders=(stk.C(1), ),
                    deleters=(stk.Br(0), ),
                ),
                stk.Bromo(
                    bromine=stk.Br(3),
                    atom=stk.C(2),
                    bonders=(stk.C(2), ),
                    deleters=(stk.Br(3), ),
                ),
            ),
            core_atoms=(stk.C(1), stk.C(2)),
            placers=(stk.C(1), stk.C(2)),
        ),
        CaseData(
            building_block=stk.BuildingBlock(
                smiles='Br[C+2][C+2]Br',
                placer_ids=(1, 2),
            ),
            functional_groups=(),
            core_atoms=(stk.Br(0), stk.C(1), stk.C(2), stk.Br(3)),
            placers=(stk.C(1), stk.C(2)),
        ),
        CaseData(
            building_block=stk.BuildingBlock(
                smiles='Br[C+2][C+2]Br',
                functional_groups=[stk.BromoFactory()],
                placer_ids=(0, 3),
            ),
            functional_groups=(
                stk.Bromo(
                    bromine=stk.Br(0),
                    atom=stk.C(1),
                    bonders=(stk.C(1), ),
                    deleters=(stk.Br(0), ),
                ),
                stk.Bromo(
                    bromine=stk.Br(3),
                    atom=stk.C(2),
                    bonders=(stk.C(2), ),
                    deleters=(stk.Br(3), ),
                ),
            ),
            core_atoms=(stk.C(1), stk.C(2)),
            placers=(stk.Br(0), stk.Br(3)),
        ),
        CaseData(
            building_block=stk.BuildingBlock(
                smiles='Br[C+2][C+2]Br',
                functional_groups=[stk.IodoFactory()],
            ),
            functional_groups=(),
            core_atoms=(stk.Br(0), stk.C(1), stk.C(2), stk.Br(3)),
            placers=(stk.Br(0), stk.C(1), stk.C(2), stk.Br(3)),
        ),
        CaseData(
            building_block=stk.BuildingBlock.init_from_file(

            ),
        ),
    ),
)
def default_init(request):
    return request.param
