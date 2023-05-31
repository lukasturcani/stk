import pytest
import stk

from ..case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda: CaseData(
            building_block=stk.BuildingBlock("Br[C+2][C+2]Br"),
            functional_groups=(),
            core_atom_ids=(0, 1, 2, 3),
            placer_ids=(0, 1, 2, 3),
        ),
        lambda: CaseData(
            building_block=stk.BuildingBlock(
                smiles="Br[C+2][C+2]Br",
                functional_groups=[stk.BromoFactory()],
            ),
            functional_groups=(
                stk.Bromo(
                    bromine=stk.Br(0),
                    atom=stk.C(1, 2),
                    bonders=(stk.C(1, 2),),
                    deleters=(stk.Br(0),),
                ),
                stk.Bromo(
                    bromine=stk.Br(3),
                    atom=stk.C(2, 2),
                    bonders=(stk.C(2, 2),),
                    deleters=(stk.Br(3),),
                ),
            ),
            core_atom_ids=(1, 2),
            placer_ids=(1, 2),
        ),
        lambda: CaseData(
            building_block=stk.BuildingBlock(
                smiles="Br[C+2][C+2]Br",
                placer_ids=(1, 2),
            ),
            functional_groups=(),
            core_atom_ids=(0, 1, 2, 3),
            placer_ids=(1, 2),
        ),
        lambda: CaseData(
            building_block=stk.BuildingBlock(
                smiles="Br[C+2][C+2]Br",
                functional_groups=[stk.BromoFactory()],
                placer_ids=(0, 3),
            ),
            functional_groups=(
                stk.Bromo(
                    bromine=stk.Br(0),
                    atom=stk.C(1, 2),
                    bonders=(stk.C(1, 2),),
                    deleters=(stk.Br(0),),
                ),
                stk.Bromo(
                    bromine=stk.Br(3),
                    atom=stk.C(2, 2),
                    bonders=(stk.C(2, 2),),
                    deleters=(stk.Br(3),),
                ),
            ),
            core_atom_ids=(1, 2),
            placer_ids=(0, 3),
        ),
        lambda: CaseData(
            building_block=stk.BuildingBlock(
                smiles="Br[C+2][C+2]Br",
                functional_groups=[stk.IodoFactory()],
            ),
            functional_groups=(),
            core_atom_ids=(0, 1, 2, 3),
            placer_ids=(0, 1, 2, 3),
        ),
        lambda: CaseData(
            building_block=stk.BuildingBlock(
                smiles="Br[C+2]Br",
                functional_groups=[stk.BromoFactory()],
            ),
            functional_groups=(
                stk.Bromo(
                    bromine=stk.Br(0),
                    atom=stk.C(1, 2),
                    bonders=(stk.C(1, 2),),
                    deleters=(stk.Br(0),),
                ),
                stk.Bromo(
                    bromine=stk.Br(2),
                    atom=stk.C(1, 2),
                    bonders=(stk.C(1, 2),),
                    deleters=(stk.Br(2),),
                ),
            ),
            core_atom_ids=(1,),
            placer_ids=(1,),
        ),
    ),
)
def default_init(request) -> CaseData:
    return request.param()
