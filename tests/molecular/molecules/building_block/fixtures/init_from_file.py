import pytest
import stk
from collections.abc import Sequence
from ..case_data import CaseData
from dataclasses import dataclass


@pytest.fixture(
    params=["building_block.mol"],
)
def path(tmpdir, request):
    return tmpdir / request.param


@dataclass(frozen=True, slots=True)
class InitFromFileData:
    building_block: stk.BuildingBlock
    init_functional_groups: Sequence[stk.FunctionalGroup]
    init_placer_ids: Sequence[int] | None
    case_data_functional_groups: Sequence[stk.FunctionalGroup]
    case_data_core_atom_ids: Sequence[int]
    case_data_placer_ids: Sequence[int]
    known_repr: str


@pytest.fixture(
    scope="session",
    params=(
        lambda: InitFromFileData(
            building_block=stk.BuildingBlock("Br[C+2][C+2]Br"),
            init_functional_groups=(),
            init_placer_ids=None,
            case_data_functional_groups=(),
            case_data_core_atom_ids=(0, 1, 2, 3),
            case_data_placer_ids=(0, 1, 2, 3),
            known_repr="BuildingBlock('Br[C+2][C+2]Br', ())",
        ),
        lambda: InitFromFileData(
            building_block=stk.BuildingBlock("Br[C+2][C+2]Br"),
            init_functional_groups=[stk.BromoFactory()],
            init_placer_ids=None,
            case_data_functional_groups=(
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
            case_data_core_atom_ids=(1, 2),
            case_data_placer_ids=(1, 2),
            known_repr=(
                "BuildingBlock('Br[C+2][C+2]Br', (Bromo(Br(0), C(1, charge=2),"
                " bonders=(C(1, charge=2),), deleters=(Br(0),)), Bromo(Br(3), "
                "C(2, charge=2), bonders=(C(2, charge=2),), deleters=(Br(3),))"
                "))"
            ),
        ),
        lambda: InitFromFileData(
            building_block=stk.BuildingBlock("Br[C+2][C+2]Br"),
            init_functional_groups=(),
            init_placer_ids=(1, 2),
            case_data_functional_groups=(),
            case_data_core_atom_ids=(0, 1, 2, 3),
            case_data_placer_ids=(1, 2),
            known_repr="BuildingBlock('Br[C+2][C+2]Br', ())",
        ),
        lambda: InitFromFileData(
            building_block=stk.BuildingBlock("Br[C+2][C+2]Br"),
            init_functional_groups=[stk.BromoFactory()],
            init_placer_ids=(0, 3),
            case_data_functional_groups=(
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
            case_data_core_atom_ids=(1, 2),
            case_data_placer_ids=(0, 3),
            known_repr=(
                "BuildingBlock('Br[C+2][C+2]Br', (Bromo(Br(0), C(1, charge=2),"
                " bonders=(C(1, charge=2),), deleters=(Br(0),)), Bromo(Br(3), "
                "C(2, charge=2), bonders=(C(2, charge=2),), deleters=(Br(3),))"
                "))"
            ),
        ),
        lambda: InitFromFileData(
            building_block=stk.BuildingBlock("Br[C+2][C+2]Br"),
            init_functional_groups=[stk.IodoFactory()],
            init_placer_ids=None,
            case_data_functional_groups=(),
            case_data_core_atom_ids=(0, 1, 2, 3),
            case_data_placer_ids=(0, 1, 2, 3),
            known_repr="BuildingBlock('Br[C+2][C+2]Br', ())",
        ),
    ),
)
def init_from_file_data(request) -> InitFromFileData:
    """
    A :class:`.InitFromFileData` instance.

    """

    return request.param()


@pytest.fixture
def init_from_file(path, init_from_file_data):
    """
    A :class:`.CaseData` instance.

    """

    data = init_from_file_data
    data.building_block.write(path)
    return CaseData(
        building_block=stk.BuildingBlock.init_from_file(
            path=path,
            functional_groups=data.init_functional_groups,
            placer_ids=data.init_placer_ids,
        ),
        functional_groups=data.case_data_functional_groups,
        core_atom_ids=data.case_data_core_atom_ids,
        placer_ids=data.case_data_placer_ids,
        known_repr=data.known_repr,
    )
