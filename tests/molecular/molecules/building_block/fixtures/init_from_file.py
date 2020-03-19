import pytest
import stk

from ..case_data import CaseData


@pytest.fixture(
    params=[
        'building_block.mol',
        'building_block.pdb',
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
        The building block, which will written to a file, so that it
        can be initialized from it.

    init_functional_groups : :class:`iterable`
        Passed to the `functional_groups` parameter of
        :meth:`.BuildingBlock.init_from_file`.

    init_placer_ids : :class:`tuple` or :class:`NoneType`
        Pass to the `placer_ids` parameter of
        :meth:`.BuildingBlock.init_from_file`.

    case_data_functional_groups : :class:`tuple`
        The functional groups the initialized building block should
        have.

    case_data_core_atoms : :class:`tuple` of :class:`.Atom`
        The core atom the initialized building block should have.

    case_data_placers : :class:`tuple` of :class:`.Atom`
        The *placer* atoms the initialized building block should
        have.

    """

    def __init__(
        self,
        building_block,
        init_functional_groups,
        init_placer_ids,
        case_data_functional_groups,
        case_data_core_atoms,
        case_data_placers,
    ):
        self.building_block = building_block
        self.init_functional_groups = init_functional_groups
        self.init_placer_ids = init_placer_ids
        self.case_data_functional_groups = case_data_functional_groups
        self.case_data_core_atoms = case_data_core_atoms
        self.case_data_placers = case_data_placers


@pytest.fixture(
    params=(
        InitFromFileData(
            building_block=stk.BuildingBlock('Br[C+2][C+2]Br'),
            init_functional_groups=(),
            init_placer_ids=None,
            case_data_functional_groups=(),
            case_data_core_atoms=(
                stk.Br(0), stk.C(1, 2), stk.C(2, 2), stk.Br(3)
            ),
            case_data_placers=(
                stk.Br(0), stk.C(1, 2), stk.C(2, 2), stk.Br(3)
            ),
        ),
    ),
)
def init_from_file_data(request):
    """
    A :class:`.InitFromFileData` instance.

    """

    return request.param


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
        core_atoms=data.case_data_core_atoms,
        placers=data.case_data_placers,
    )
