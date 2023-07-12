import pytest
import stk

from ..case_data import CaseData


@pytest.fixture(
    params=[
        "building_block.mol",
    ],
)
def path(tmpdir, request):
    return str(tmpdir / request.param)


class InitFromFileData:
    """
    Stores data for the :meth:`.BuildingBlock.init_from_file`.

    Attributes
    ----------
    building_block : :class:`.BuildingBlock`
        The building block, which will be written to a file, so that it
        can be initialized from it.

    init_functional_groups : :class:`iterable`
        Passed to the `functional_groups` parameter of
        :meth:`.BuildingBlock.init_from_file`.

    init_placer_ids : :class:`tuple` or :class:`NoneType`
        Passed to the `placer_ids` parameter of
        :meth:`.BuildingBlock.init_from_file`.

    case_data_functional_groups : :class:`tuple`
        The functional groups the initialized building block should
        have.

    case_data_core_atom_ids : :class:`tuple` of :class:`int`
        The ids of core atoms the initialized building block should
        have.

    case_data_placer_ids : :class:`tuple` of :class:`int`
        The ids of *placer* atoms the initialized building block should
        have.

    """

    def __init__(
        self,
        building_block,
        init_functional_groups,
        init_placer_ids,
        case_data_functional_groups,
        case_data_core_atom_ids,
        case_data_placer_ids,
    ):
        """
        Initialize a :class:`.InitFromData` instance.

        Parameters
        ----------
        building_block : :class:`.BuildingBlock`
            The building block, which will be written to a file, so
            that it can be initialized from it.

        init_functional_groups : :class:`iterable`
            Passed to the `functional_groups` parameter of
            :meth:`.BuildingBlock.init_from_file`.

        init_placer_ids : :class:`tuple` or :class:`NoneType`
            Passed to the `placer_ids` parameter of
            :meth:`.BuildingBlock.init_from_file`.

        case_data_functional_groups : :class:`tuple`
            The functional groups the initialized building block should
            have.

        case_data_core_atom_ids : :class:`tuple` of :class:`int`
            The ids of core atoms the initialized building block should
            have.

        case_data_placer_ids : :class:`tuple` of :class:`int`
            The ids of *placer* atoms the initialized building block
            should have.

        """

        self.building_block = building_block
        self.init_functional_groups = init_functional_groups
        self.init_placer_ids = init_placer_ids
        self.case_data_functional_groups = case_data_functional_groups
        self.case_data_core_atom_ids = case_data_core_atom_ids
        self.case_data_placer_ids = case_data_placer_ids


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
        ),
        lambda: InitFromFileData(
            building_block=stk.BuildingBlock("Br[C+2][C+2]Br"),
            init_functional_groups=(),
            init_placer_ids=(1, 2),
            case_data_functional_groups=(),
            case_data_core_atom_ids=(0, 1, 2, 3),
            case_data_placer_ids=(1, 2),
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
        ),
        lambda: InitFromFileData(
            building_block=stk.BuildingBlock("Br[C+2][C+2]Br"),
            init_functional_groups=[stk.IodoFactory()],
            init_placer_ids=None,
            case_data_functional_groups=(),
            case_data_core_atom_ids=(0, 1, 2, 3),
            case_data_placer_ids=(0, 1, 2, 3),
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
    )
