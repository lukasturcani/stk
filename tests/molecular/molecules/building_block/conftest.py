import stk
import itertools as it
import pytest


from .case_data import CaseData


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
        self.case_functional_groups = case_data_functional_groups
        self.case_core_atoms = case_data_core_atoms
        self.case_placers = case_data_placers


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
    return request.param


@pytest.fixture
def case_data_2(path, init_from_file_data):
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
def case_data(request):
    """
    A :class:`.CaseData` instance.

    """

    return request.param.clone()


@pytest.fixture(
    params=(
        lambda molecule:
            stk.BromoFactory().get_functional_groups(molecule),
        lambda molecule:
            stk.PrimaryAminoFactory().get_functional_groups(molecule),
        lambda molecule: it.chain(
            stk.PrimaryAminoFactory().get_functional_groups(molecule),
            stk.BromoFactory().get_functional_groups(molecule)),
    )
)
def get_functional_groups(request):
    """
    Yield the functional groups of a `molecule`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule whose functional groups should be gotten.

    Yields
    ------
    :class:`.FunctionalGroup`
        A functional group of `molecule`.

    """

    return request.param
