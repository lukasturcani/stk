import pytest
import stk

from ..case_data import CaseData


class CofData:
    """
    Data to initialize a COF :class:`.CaseData` instance.

    Attributes
    ----------
    topology_graph : :class:`type`
        A COF class.

    building_blocks : :class:`tuple` of :class:`.BuildingBlock`
        The building blocks of the COF.

    lattice_size : :class:`tuple` of :class:`int`
        The size of the lattice.

    vertex_alignments : :class:`dict`
        Passed to the `vertex_alignments` parameter of the COF
        initializer.

    num_new_atoms : :class:`int`
        The number of new atoms added by the construction process.

    num_new_bonds : :class:`int`
        The number of new bonds added by the construction process.

    num_building_blocks : :class:`dict`
        For each building block in :attr:`building_blocks`, maps its
        index to the number of times its used in the construction of
        the COF.

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
        Initialize a :class:`.CofData` instance.

        Parameters
        ----------
        topology_graph : :class:`type`
            A COF class.

        building_blocks : :class:`tuple` of :class:`.BuildingBlock`
            The building blocks of the COF.

        lattice_size : :class:`tuple` of :class:`int`
            The size of the lattice.

        vertex_alignments : :class:`dict`
            Passed to the `vertex_alignments` parameter of the COF
            initializer.

        num_new_atoms : :class:`int`
            The number of new atoms added by the construction process.

        num_new_bonds : :class:`int`
            The number of new bonds added by the construction process.

        num_building_blocks : :class:`dict`
            For each building block in `building_blocks`, maps its
            index to the number of times its used in the construction
            of the COF.

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

    @classmethod
    def init_from_construction_result(
        cls,
        topology_graph,
        building_blocks,
        lattice_size,
        vertex_alignments,
        num_new_atoms,
        num_new_bonds,
        num_building_blocks,
    ):
        """
        Initialize a :class:`.CofData` instance.

        This method creates the constructed molecule using
        :meth:`.ConstructedMolecule.init_from_construction_result`.

        Parameters
        ----------
        topology_graph : :class:`type`
            A COF class.

        building_blocks : :class:`tuple` of :class:`.BuildingBlock`
            The building blocks of the COF.

        lattice_size : :class:`tuple` of :class:`int`
            The size of the lattice.

        vertex_alignments : :class:`dict`
            Passed to the `vertex_alignments` parameter of the COF
            initializer.

        num_new_atoms : :class:`int`
            The number of new atoms added by the construction process.

        num_new_bonds : :class:`int`
            The number of new bonds added by the construction process.

        num_building_blocks : :class:`dict`
            For each building block in `building_blocks`, maps its
            index to the number of times its used in the construction
            of the COF.

        """

        obj = cls.__new__(cls)
        topology_graph_instance = topology_graph(
            building_blocks=building_blocks,
            lattice_size=lattice_size,
            vertex_alignments=vertex_alignments,
        )
        construction_result = topology_graph_instance.construct()
        obj.constructed_molecule = (
            stk.ConstructedMolecule.init_from_construction_result(
                construction_result=construction_result,
            )
        )
        obj.num_new_atoms = num_new_atoms
        obj.num_new_bonds = num_new_bonds
        obj.num_building_blocks = {
            building_blocks[index]: num
            for index, num in num_building_blocks.items()
        }
        obj.building_blocks = building_blocks
        return obj


@pytest.fixture(
    scope="session",
    params=(
        lambda: CofData(
            topology_graph=stk.cof.Honeycomb,
            building_blocks=(
                stk.BuildingBlock(
                    smiles=("Br[C+]1[C+2][C+](Br)[C+](F)[C+](Br)[C+2]1"),
                    functional_groups=[stk.BromoFactory()],
                ),
                stk.BuildingBlock(
                    smiles="Br[C+]=NC#CBr",
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
            lattice_size=(2, 2, 1),
            vertex_alignments=None,
            num_new_atoms=0,
            num_new_bonds=20,
            num_building_blocks={0: 8, 1: 12},
        ),
        lambda: CofData.init_from_construction_result(
            topology_graph=stk.cof.Honeycomb,
            building_blocks=(
                stk.BuildingBlock(
                    smiles=("Br[C+]1[C+2][C+](Br)[C+](F)[C+](Br)[C+2]1"),
                    functional_groups=[stk.BromoFactory()],
                ),
                stk.BuildingBlock(
                    smiles="Br[C+]=NC#CBr",
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
            lattice_size=(2, 2, 1),
            vertex_alignments=None,
            num_new_atoms=0,
            num_new_bonds=20,
            num_building_blocks={0: 8, 1: 12},
        ),
        lambda: CofData(
            topology_graph=stk.cof.Honeycomb,
            building_blocks=(
                stk.BuildingBlock(
                    smiles=("Br[C+]1[C+2][C+](Br)[C+](F)[C+](Br)[C+2]1"),
                    functional_groups=[stk.BromoFactory()],
                ),
                stk.BuildingBlock(
                    smiles="Br[C+]=NC#CBr",
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
            lattice_size=(2, 2, 1),
            vertex_alignments={0: 1, 1: 1, 2: 1, 3: 1, 4: 1},
            num_new_atoms=0,
            num_new_bonds=20,
            num_building_blocks={0: 8, 1: 12},
        ),
        lambda: CofData(
            topology_graph=stk.cof.Honeycomb,
            building_blocks=(
                stk.BuildingBlock(
                    smiles=("Br[C+]1[C+2][C+](Br)[C+](F)[C+](Br)[C+2]1"),
                    functional_groups=[stk.BromoFactory()],
                ),
                stk.BuildingBlock(
                    smiles="Br[C+]=NC#CBr",
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
            lattice_size=(2, 2, 1),
            vertex_alignments={0: 2, 1: 2},
            num_new_atoms=0,
            num_new_bonds=20,
            num_building_blocks={0: 8, 1: 12},
        ),
        lambda: CofData(
            topology_graph=stk.cof.LinkerlessHoneycomb,
            building_blocks=(
                stk.BuildingBlock(
                    smiles=("Br[C+]1[C+2][C+](Br)[C+](F)[C+](Br)[C+2]1"),
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
            lattice_size=(2, 2, 1),
            vertex_alignments=None,
            num_new_atoms=0,
            num_new_bonds=8,
            num_building_blocks={0: 8},
        ),
        lambda: CofData(
            topology_graph=stk.cof.Hexagonal,
            building_blocks=(
                stk.BuildingBlock(
                    smiles=("Br[C+]1[C+](F)[C+](I)[C+](I)[C+](Br)C1Br"),
                    functional_groups=[
                        stk.BromoFactory(),
                        stk.IodoFactory(),
                        stk.FluoroFactory(),
                    ],
                ),
                stk.BuildingBlock(
                    smiles="Br[C+]=NC#CBr",
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
            lattice_size=(2, 2, 1),
            vertex_alignments={0: 5},
            num_new_atoms=0,
            num_new_bonds=81,
            num_building_blocks={0: 16, 1: 48},
        ),
        lambda: CofData(
            topology_graph=stk.cof.Kagome,
            building_blocks=(
                stk.BuildingBlock(
                    smiles=("Br[C+]1[C+](Br)[C+](F)[C+](Br)[C+](Br)" "[C+2]1"),
                    functional_groups=[stk.BromoFactory()],
                ),
                stk.BuildingBlock(
                    smiles="Br[C+]=NC#CBr",
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
            lattice_size=(2, 2, 1),
            vertex_alignments=None,
            num_new_atoms=0,
            num_new_bonds=41,
            num_building_blocks={0: 12, 1: 24},
        ),
        lambda: CofData(
            topology_graph=stk.cof.Square,
            building_blocks=(
                stk.BuildingBlock(
                    smiles="BrC1=C(Br)C(F)(Br)[C+]1Br",
                    functional_groups=[stk.BromoFactory()],
                ),
                stk.BuildingBlock(
                    smiles="Br[C+]=NC#CBr",
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
def cof_data(request) -> CofData:
    """
    A :class:`.CofData` instance.

    """

    return request.param()


@pytest.fixture
def cof(cof_data):
    """
    A :class:`.CaseData` instance.

    """

    return CaseData(
        constructed_molecule=cof_data.constructed_molecule,
        num_new_atoms=cof_data.num_new_atoms,
        num_new_bonds=cof_data.num_new_bonds,
        num_building_blocks=cof_data.num_building_blocks,
        building_blocks=cof_data.building_blocks,
    )
