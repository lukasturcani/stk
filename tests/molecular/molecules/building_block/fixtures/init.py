import numpy as np
import pytest
import stk

from ..case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda: CaseData(
            building_block=stk.BuildingBlock.init(
                atoms=(stk.Br(0), stk.C(1, 2), stk.C(2, 2), stk.Br(3)),
                bonds=(
                    stk.Bond(stk.Br(0), stk.C(1, 2), 1),
                    stk.Bond(stk.C(1, 2), stk.C(2, 2), 1),
                    stk.Bond(stk.C(2, 2), stk.Br(3), 1),
                ),
                position_matrix=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [1.0, 1.0, 1.0],
                        [2.0, 2.0, 2.0],
                        [3.0, 3.0, 3.0],
                    ]
                ),
            ),
            functional_groups=(),
            core_atom_ids=(0, 1, 2, 3),
            placer_ids=(0, 1, 2, 3),
        ),
        lambda: CaseData(
            building_block=stk.BuildingBlock.init(
                atoms=(stk.Br(0), stk.C(1, 2), stk.C(2, 2), stk.Br(3)),
                bonds=(
                    stk.Bond(stk.Br(0), stk.C(1, 2), 1),
                    stk.Bond(stk.C(1, 2), stk.C(2, 2), 1),
                    stk.Bond(stk.C(2, 2), stk.Br(3), 1),
                ),
                position_matrix=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [1.0, 1.0, 1.0],
                        [2.0, 2.0, 2.0],
                        [3.0, 3.0, 3.0],
                    ]
                ),
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
            building_block=stk.BuildingBlock.init(
                atoms=(stk.Br(0), stk.C(1, 2), stk.C(2, 2), stk.Br(3)),
                bonds=(
                    stk.Bond(stk.Br(0), stk.C(1, 2), 1),
                    stk.Bond(stk.C(1, 2), stk.C(2, 2), 1),
                    stk.Bond(stk.C(2, 2), stk.Br(3), 1),
                ),
                position_matrix=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [1.0, 1.0, 1.0],
                        [2.0, 2.0, 2.0],
                        [3.0, 3.0, 3.0],
                    ]
                ),
                placer_ids=(1, 2),
            ),
            functional_groups=(),
            core_atom_ids=(0, 1, 2, 3),
            placer_ids=(1, 2),
        ),
        lambda: CaseData(
            building_block=stk.BuildingBlock.init(
                atoms=(stk.Br(0), stk.C(1, 2), stk.C(2, 2), stk.Br(3)),
                bonds=(
                    stk.Bond(stk.Br(0), stk.C(1, 2), 1),
                    stk.Bond(stk.C(1, 2), stk.C(2, 2), 1),
                    stk.Bond(stk.C(2, 2), stk.Br(3), 1),
                ),
                position_matrix=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [1.0, 1.0, 1.0],
                        [2.0, 2.0, 2.0],
                        [3.0, 3.0, 3.0],
                    ]
                ),
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
            building_block=stk.BuildingBlock.init(
                atoms=(stk.Br(0), stk.C(1, 2), stk.C(2, 2), stk.Br(3)),
                bonds=(
                    stk.Bond(stk.Br(0), stk.C(1, 2), 1),
                    stk.Bond(stk.C(1, 2), stk.C(2, 2), 1),
                    stk.Bond(stk.C(2, 2), stk.Br(3), 1),
                ),
                position_matrix=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [1.0, 1.0, 1.0],
                        [2.0, 2.0, 2.0],
                        [3.0, 3.0, 3.0],
                    ]
                ),
                functional_groups=[stk.IodoFactory()],
            ),
            functional_groups=(),
            core_atom_ids=(0, 1, 2, 3),
            placer_ids=(0, 1, 2, 3),
        ),
    ),
)
def init(request) -> CaseData:
    return request.param()
