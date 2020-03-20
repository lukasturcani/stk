import pytest
import stk

from ..case_data import CaseData


class LinearData:
    """
    Data to initialize a linear polymer :class:`.CaseData` instance.

    Attributes
    ----------
    building_blocks : :class:`tuple` of :class:`.BuildingBlock`
        The building blocks of the polymer.

    repeating_unit : :class:`str` or :class:`tuple` of :class:`int`
        Passed to the `repeating_unit` parameter of
        :meth:`.polymer.Linear.__init__`.

    num_repeating_units : :class:`int`
        Passed to the `num_repeating_units` parameter of
        :meth:`.polymer.Linear.__init__`.

    num_new_atoms : :class:`int`
        The number of new atoms added by the construction process.

    num_new_bonds : :class:`int`
        The number of new bonds added by the construction process.

    num_building_blocks : :class:`dict`
        For each building block in :attr:`building_blocks`, maps its
        index to the number of times its used in the construction of
        the polymer.

    """

    def __init__(
        self,
        building_blocks,
        repeating_unit,
        num_repeating_units,
        num_new_atoms,
        num_new_bonds,
        num_building_blocks,
    ):
        """
        Initiailze a :class:`.LinearData` instance.

        Parameters
        ----------
        building_blocks : :class:`tuple` of :class:`.BuildingBlock`
            The building blocks of the polymer.

        repeating_unit : :class:`str` or :class:`tuple` of :class:`int`
            Passed to the `repeating_unit` parameter of
            :meth:`.polymer.Linear.__init__`.

        num_repeating_units : :class:`int`
            Passed to the `num_repeating_units` parameter of
            :meth:`.polymer.Linear.__init__`.

        num_new_atoms : :class:`int`
            The number of new atoms added by the construction process.

        num_new_bonds : :class:`int`
            The number of new bonds added by the construction process.

        num_building_blocks : :class:`dict`
            For each building block in `building_blocks`, maps its
            index to the number of times its used in the construction
            of the polymer.

        """

        self.building_blocks = building_blocks
        self.repeating_unit = repeating_unit
        self.num_repeating_units = num_repeating_units
        self.num_new_atoms = num_new_atoms
        self.num_new_bonds = num_new_bonds
        self.num_building_blocks = num_building_blocks


@pytest.fixture(
    params=(
        LinearData(
            building_blocks=(
                stk.BuildingBlock('BrC#CBr', [stk.BromoFactory()]),
                stk.BuildingBlock(
                    smiles='Br[C+]=NC#CBr',
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
            repeating_unit='AB',
            num_repeating_units=2,
            num_new_atoms=0,
            num_new_bonds=3,
        ),
    ),
)
def linear_data(request):
    """
    A :class:`.LinearData` instance.

    """

    return request.param


@pytest.fixture
def linear(linear_data):
    """
    A :class:`.CaseData` instance.

    """

    return CaseData(
        topology_graph=stk.polymer.Linear(
            building_blocks=linear_data.building_blocks,
            repeating_unit=linear_data.repeating_unit,
            num_repeating_units=linear_data.num_repeating_units,
        ),
        num_new_atoms=linear_data.num_new_atoms,
        num_new_bonds=linear_data.num_new_bonds,
        num_building_blocks={
            linear_data.building_blocks[index]: num
            for index, num in linear_data.num_building_blocks.items()
        },
        building_blocks=linear_data.building_blocks,
    )
