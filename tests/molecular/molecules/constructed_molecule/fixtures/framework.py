import pytest
import stk

from ..case_data import CaseData


class FrameworkData:
    """
    Data to initialize a framework :class:`.CaseData` instance.

    Attributes
    ----------
    topology_graph : :class:`type`
        A framework class.

    building_blocks : :class:`tuple` of :class:`.BuildingBlock`
        The building blocks of the framework.

    lattice_size : :class:`tuple` of :class:`int`
        The size of the lattice.

    vertex_alignments : :class:`dict`
        Passed to the `vertex_alignments` parameter of the framework
        initializer.

    num_new_atoms : :class:`int`
        The number of new atoms added by the construction process.

    num_new_bonds : :class:`int`
        The number of new bonds added by the construction process.

    num_building_blocks : :class:`dict`
        For each building block in :attr:`building_blocks`, maps its
        index to the number of times its used in the construction of
        the framework.

    """

    def __init__(
        self,
        topology_graph,
        building_blocks,
        lattice_size,
        vertex_alignments,
        num_new_atoms,
        num_new_bonds,
        num_building_blocks,
    ):
        """
        Initialize a :class:`.FrameworkData` instance.

        Parameters
        ----------
        topology_graph : :class:`type`
            A framework class.

        building_blocks : :class:`tuple` of :class:`.BuildingBlock`
            The building blocks of the framework.

        lattice_size : :class:`tuple` of :class:`int`
            The size of the lattice.

        vertex_alignments : :class:`dict`
            Passed to the `vertex_alignments` parameter of the
            framework initializer.

        num_new_atoms : :class:`int`
            The number of new atoms added by the construction process.

        num_new_bonds : :class:`int`
            The number of new bonds added by the construction process.

        num_building_blocks : :class:`dict`
            For each building block in `building_blocks`, maps its
            index to the number of times its used in the construction
            of the framework.

        """

        self.constructed_molecule = stk.ConstructedMolecule(
            topology_graph=topology_graph(
                building_blocks=building_blocks,
                lattice_size=lattice_size,
                vertex_alignments=vertex_alignments,
            )
        )
        self.num_new_atoms = num_new_atoms
        self.num_new_bonds = num_new_bonds
        self.num_building_blocks = {
            building_blocks[index]: num
            for index, num in num_building_blocks.items()
        }
        self.building_blocks = building_blocks


@pytest.fixture(
    params=(
        FrameworkData(
            topology_graph=stk.framework.Honeycomb,
            building_blocks=(
                stk.BuildingBlock(
                    smiles=(
                        'Br[C+]1[C+2][C+](Br)[C+](F)[C+](Br)[C+2]1'
                    ),
                    functional_groups=[stk.BromoFactory()],
                ),
                stk.BuildingBlock(
                    smiles='Br[C+]=NC#CBr',
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
            lattice_size=(2, 2, 1),
            vertex_alignments=None,
            num_new_atoms=0,
            num_new_bonds=20,
            num_building_blocks={0: 8, 1: 12},
        ),
        FrameworkData(
            topology_graph=stk.framework.Honeycomb,
            building_blocks=(
                stk.BuildingBlock(
                    smiles=(
                        'Br[C+]1[C+2][C+](Br)[C+](F)[C+](Br)[C+2]1'
                    ),
                    functional_groups=[stk.BromoFactory()],
                ),
                stk.BuildingBlock(
                    smiles='Br[C+]=NC#CBr',
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
            lattice_size=(2, 2, 1),
            vertex_alignments={0: 1, 1: 1, 2: 1, 3: 1, 4: 1},
            num_new_atoms=0,
            num_new_bonds=20,
            num_building_blocks={0: 8, 1: 12},
        ),
        FrameworkData(
            topology_graph=stk.framework.Honeycomb,
            building_blocks=(
                stk.BuildingBlock(
                    smiles=(
                        'Br[C+]1[C+2][C+](Br)[C+](F)[C+](Br)[C+2]1'
                    ),
                    functional_groups=[stk.BromoFactory()],
                ),
                stk.BuildingBlock(
                    smiles='Br[C+]=NC#CBr',
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
            lattice_size=(2, 2, 1),
            vertex_alignments={0: 2, 1: 2},
            num_new_atoms=0,
            num_new_bonds=20,
            num_building_blocks={0: 8, 1: 12},
        ),
        FrameworkData(
            topology_graph=stk.framework.LinkerlessHoneycomb,
            building_blocks=(
                stk.BuildingBlock(
                    smiles=(
                        'Br[C+]1[C+2][C+](Br)[C+](F)[C+](Br)[C+2]1'
                    ),
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
            lattice_size=(2, 2, 1),
            vertex_alignments=None,
            num_new_atoms=0,
            num_new_bonds=8,
            num_building_blocks={0: 8},
        ),
        FrameworkData(
            topology_graph=stk.framework.Hexagonal,
            building_blocks=(
                stk.BuildingBlock(
                    smiles=(
                        'Br[C+]1[C+](F)[C+](I)[C+](I)[C+](Br)C1Br'
                    ),
                    functional_groups=[
                        stk.BromoFactory(),
                        stk.IodoFactory(),
                        stk.FluoroFactory(),
                    ],
                ),
                stk.BuildingBlock(
                    smiles='Br[C+]=NC#CBr',
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
            lattice_size=(2, 2, 1),
            vertex_alignments={0: 5},
            num_new_atoms=0,
            num_new_bonds=81,
            num_building_blocks={0: 16, 1: 48},
        ),
        FrameworkData(
            topology_graph=stk.framework.Kagome,
            building_blocks=(
                stk.BuildingBlock(
                    smiles=(
                        'Br[C+]1[C+](Br)[C+](F)[C+](Br)[C+](Br)'
                        '[C+2]1'
                    ),
                    functional_groups=[stk.BromoFactory()],
                ),
                stk.BuildingBlock(
                    smiles='Br[C+]=NC#CBr',
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
            lattice_size=(2, 2, 1),
            vertex_alignments=None,
            num_new_atoms=0,
            num_new_bonds=41,
            num_building_blocks={0: 12, 1: 24},
        ),
        FrameworkData(
            topology_graph=stk.framework.Square,
            building_blocks=(
                stk.BuildingBlock(
                    smiles='BrC1=C(Br)C(F)(Br)[C+]1Br',
                    functional_groups=[stk.BromoFactory()],
                ),
                stk.BuildingBlock(
                    smiles='Br[C+]=NC#CBr',
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
            lattice_size=(2, 2, 1),
            vertex_alignments=None,
            num_new_atoms=0,
            num_new_bonds=12,
            num_building_blocks={0: 4, 1: 8},
        ),
    ),
)
def framework_data(request):
    """
    A :class:`.FrameworkData` instance.

    """

    return request.param


@pytest.fixture
def framework(framework_data):
    """
    A :class:`.CaseData` instance.

    """

    return CaseData(
        constructed_molecule=framework_data.constructed_molecule,
        num_new_atoms=framework_data.num_new_atoms,
        num_new_bonds=framework_data.num_new_bonds,
        num_building_blocks=framework_data.num_building_blocks,
        building_blocks=framework_data.building_blocks,
    )
