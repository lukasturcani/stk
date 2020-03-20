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

        self.constructed_molecule = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=building_blocks,
                repeating_unit=repeating_unit,
                num_repeating_units=num_repeating_units,
            ),
        )
        self.num_new_atoms = num_new_atoms
        self.num_new_bonds = num_new_bonds
        self.num_building_blocks = {
            building_blocks[index]: num
            for index, num in num_building_blocks.items()
        }
        self.building_blocks = building_blocks

    @classmethod
    def init_with_building_blocks(cls, building_blocks):
        """
        Initialize by substituing building blocks.

        Initializes a case by running
        :meth:`.ConstructedMolecule.with_building_blocks`.

        Parameters
        ----------
        building_blocks : :class:`tuple` of :class:`.BuildingBlock`
            The building blocks the polymer being tested should hold.

        """

        initial = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrCNC#CBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles='Br[C+]=NCNC#CBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                repeating_unit='AB',
                num_repeating_units=3,
            )
        )
        building_block_map = {
            building_block1: building_block2
            for building_block1, building_block2
            in zip(
                initial.get_building_blocks(),
                building_blocks,
            )
        }

        obj = cls.__new__(cls)
        obj.constructed_molecule = initial.with_building_blocks(
            building_block_map=building_block_map,
        )
        obj.num_new_atoms = 0
        obj.num_new_bonds = 5
        obj.building_blocks = building_blocks
        obj.num_building_blocks = {
            building_block2: initial.get_num_building_block(
                building_block=building_block1,
            )
            for building_block1, building_block2
            in zip(
                initial.get_building_blocks(),
                building_blocks,
            )
        }
        return obj


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
            num_building_blocks={0: 2, 1: 2},
        ),
        LinearData.init_with_building_blocks(
            building_blocks=(
                stk.BuildingBlock('BrC#CBr', [stk.BromoFactory()]),
                stk.BuildingBlock(
                    smiles='Br[C+]=NC#CBr',
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
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
        constructed_molecule=linear_data.constructed_molecule,
        num_new_atoms=linear_data.num_new_atoms,
        num_new_bonds=linear_data.num_new_bonds,
        num_building_blocks=linear_data.num_building_blocks,
        building_blocks=linear_data.building_blocks,
    )
