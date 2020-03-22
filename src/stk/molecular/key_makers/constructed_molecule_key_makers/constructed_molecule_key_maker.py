"""
Constructed Molecule Key Maker
==============================

"""

from .topology_graph_key_makers import (
    CageKeyMaker,
    CofKeyMaker,
    LinearPolymerKeyMaker,
    HostGuestComplexKeyMaker,
    NRotaxaneKeyMaker,
    MacrocycleKeyMaker,
)
from ..molecule_key_makers import InchiKey
from ...topology_graphs.cage import (
    EightPlusSixteen,
    EightPlusTwelve,
    FivePlusTen,
    FourPlusEight,
    FourPlusFour,
    FourPlusSix,
    FourPlusSix2,
    OnePlusOne,
    SixPlusEight,
    SixPlusNine,
    SixPlusTwelve,
    TenPlusTwenty,
    ThreePlusSix,
    TwelvePlusThirty,
    TwentyPlusThirty,
    TwoPlusFour,
    TwoPlusThree,
    TwoPlusTwo,
)
from ..topology_graphs.cof import (
    Hexagonal,
    Honeycomb,
    Kagome,
    LinkerlessHoneycomb,
    Square,
)
from ..topology_graphs.polymer import Linear
from ..topology_graphs.host_guest import Complex
from ..topology_graphs.rotaxane import NRotaxane
from ..topology_graphs.macrocycle import Macrocycle


class ConstructedMoleculeKeyMaker:
    """
    Abstract base class for making :class:`.ConstructedMolecule` keys.

    Notes
    -----
    You might notice that the public methods of this abstract base
    class are implemented. This is purely for the convenience of users.
    The implemented public methods are simply default implementations,
    which can be safely ignored or overridden, when implementing
    subclasses. However, the default implementation works without
    need for any further subclassing, and can be used directly, if it
    suits your needs.

    Examples
    --------
    *Subclass Implementation*

    There's just two simple methods to implement.

    .. code-block:: python

        import stk


        class NumBuildingBlocks(stk.ConstructedMoleculeKeyMaker):
            def __init__(self):
                return

            def name(self):
                # What string this is is completely up to you. It does
                # not have to be related to the class name.
                return 'num_building_blocks'

            def get_key(self, constructed_molecule):
                building_blocks = (
                    constructed_molecule.get_building_blocks()
                )
                return sum(1 for _ in building_blocks)


        # A usage example of the new subclass.

        jsonizer = stk.ConstructedMoleculeJsonizer(
            constructed_molecule_key_makers=(NumBuildingBlocks(), ),
        )

        # Make a constructed molecule, which you want to convert to a
        # JSON format.

        building_block = stk.BuildingBlock(
            smiles='BrCCBr',
            functional_groups=[stk.BromoFactory()],
        )
        polymer_topology_graph = stk.polymer.Linear(
            building_blocks=(building_block, ),
            repeating_unit='A',
            num_repeating_units=12,
        )
        polymer = stk.ConstructedMolecule(polymer_topology_graph)

        # Create a JSON representation of a constructed molecule, which
        # holds the number of building blocks.
        json = jsonizer.to_json(polymer)

    *Usage*

    Because :class:`.ConstructedMoleculeKeyMaker` comes with a default
    implementation, it can be used directly, instead of having to
    make a subclass

    .. code-block:: python

        import stk

        jsonizer = stk.ConstructedMoleculeJsonizer(
            constructed_molecule_key_makers=(
                stk.ConstructedMoleculeKeyMaker(),
            ),
        )
        # Create a JSON representation of a constructed molecule, which
        # holds the default constructed molecule key.
        json = jsonizer.to_json(polymer)

    You want to use an InChI instead of an InChIKey in the
    constructed molecule key

    .. code-block:: python

        jsonizer = stk.ConstructedMoleculeJsonizer(
            constructed_molecule_key_makers=(
                stk.ConstructedMoleculeKeyMaker(
                    # If you change the nature of the key, its a good
                    # idea to change its name to reflect that.
                    name='InChIConstructedMoleculeKey',
                    molecule_key_maker=stk.Inchi(),
                ),
                # You can still keep the default constructed molecule
                # key maker too. No pressure though, excluding it from
                # this tuple is also valid.
                stk.ConstructedMoleculeKeyMaker(),
            ),
        )
        # Create a JSON representation of a constructed molecule, which
        # holds two constructed molecule keys, one featuring an
        # InChI and one featuring an InChIKey.
        json = jsonizer.to_json(polymer)

    You want to use a different topology graph key maker for some
    :class:`.TopologyGraph` subclass

    .. code-block:: python

        import stk

        # Define the alternate topology graph key maker you want to
        # use. In this case, the polymer key is given by the number of
        # building blocks.
        class LinearKeyMaker:
            def get_key(self, linear):
                building_blocks = (
                    linear.get_building_blocks()
                )
                return sum(1 for _ in building_blocks)

        key_maker = stk.ConstructedMoleculeKeyMaker(
            # Use an alternate topology graph key maker for
            # stk.polymer.Linear instances.
            topology_graph_key_makers={
                stk.polymer.Linear: LinearKeyMaker(),
            }
        )

        # Use the key_maker as normal.
        jsonizer = stk.ConstructedMoleculeJsonizer(
            constructed_molecule_key_makers=(key_maker, ),
        )
        json = jsonizer.to_json(polymer)

    If you want to use a key maker for a new :class:`.TopologyGraph`
    the process is essentially the same as the one described above

    .. code-block:: python

        # Define your TopologyGraph subclass
        class MyTopologyGraph(stk.TopologyGraph):
            # Your implementation goes here


        # Define the topology graph key maker for your subclass
        class MyTopologyGraphKeyMaker:
            def get_key(self, my_topology_graph):
                # Get the key for my_topology_graph somehow.


        key_maker = stk.ConstructedMoleculeKeyMaker(
            # use the new key maker.
            topology_graph_key_makers={
                MyTopologyGraph: MyTopologyGraphKeyMaker(),
            },
        )

    """

    def __init__(
        self,
        name='ConstructedMoleculeKey',
        molecule_key_maker=InchiKey(),
        topology_graph_key_makers=None,
    ):
        """
        Initialize a :class:`.ConstructedMoleculeKeyMaker` instance.

        Parameters
        ----------
        name : :class:`str`, optional
            The name of the key made by the maker.

        molecule_key_maker : :class:`.MoleculeKeyMaker`, optional
            Used to generate the part of the key responsible for the
            molecular component of :class:`.ConstructedMolecule`.

        topology_graph_key_makers : :class:`dict`, optional
            Maps a :class:`.TopologyGraph` subclass to the
            topology graph key maker, which should be used to get key
            for its instances. For a list of built-in topology graph
            key makers, see :mod:`.topology_graph_key_makers`. If
            ``None``, the built-in key makers will be used.

            This parameter will only update the default mapping.
            Providing an empty ``{}``, will mean that all the
            defaults are kept. Providing

            .. code-block:: python

                key_maker = stk.ConstructedMoleculeKeyMaker(
                    topology_graph_key_makers={
                        stk.polymer.Linear: CustomLinearKeyMaker(),
                    },
                )

            will keep all the defaults, with the exception of
            :class:`.polymer.Linear`, which will use
            :class:`.CustomLinearKeyMaker` instead of the default.

        """

        if topology_graph_key_makers is None:
            topology_graph_key_makers = {}

        self._name = name
        self._molecule_key_maker = molecule_key_maker
        self._topology_graph_key_makers = (
            self._get_default_topology_graph_key_makers()
        )
        self._topology_graph_key_makers.update(
            topology_graph_key_makers
        )
        # Used for __repr__().
        self._input_topology_graph_key_makers = (
            topology_graph_key_makers
        )

    def get_name(self):
        """
        Get the name of the key.

        Returns
        -------
        :class:`str`
            The name.

        """

        return self._name

    def get_key(self, constructed_molecule):
        """
        Get the key of `constructed_molecule`.

        Parameters
        ----------
        constructed_molecule : :class:`.ConstructedMolecule`
            The constructed molecule for which a key is wanted.

        Returns
        -------
        :class:`object`

        """

        # Use private attribute access because accessing the topology
        # graph is not part of the public interface of a constructed
        # molecule.
        topology_graph = constructed_molecule._topology_graph
        key_maker = (
            self._topology_graph_key_makers[type(topology_graph)]
        )
        return key_maker.get_key(topology_graph)

    @staticmethod
    def _get_default_topology_graph_keys():
        """
        Get the default topology graph key makers.

        Returns
        -------
        :class:`dict`
            Maps a :class:`.TopologyGrpah` subclass to the
            topology graph key maker that should be used with it.

        """

        return {
            EightPlusSixteen: CageKeyMaker(),
            EightPlusTwelve: CageKeyMaker(),
            FivePlusTen: CageKeyMaker(),
            FourPlusEight: CageKeyMaker(),
            FourPlusFour: CageKeyMaker(),
            FourPlusSix: CageKeyMaker(),
            FourPlusSix2: CageKeyMaker(),
            OnePlusOne: CageKeyMaker(),
            SixPlusEight: CageKeyMaker(),
            SixPlusNine: CageKeyMaker(),
            SixPlusTwelve: CageKeyMaker(),
            TenPlusTwenty: CageKeyMaker(),
            ThreePlusSix: CageKeyMaker(),
            TwelvePlusThirty: CageKeyMaker(),
            TwentyPlusThirty: CageKeyMaker(),
            TwoPlusFour: CageKeyMaker(),
            TwoPlusThree: CageKeyMaker(),
            TwoPlusTwo: CageKeyMaker(),
            Hexagonal: CofKeyMaker(),
            Honeycomb: CofKeyMaker(),
            Kagome: CofKeyMaker(),
            LinkerlessHoneycomb: CofKeyMaker(),
            Square: CofKeyMaker(),
            Linear: LinearPolymerKeyMaker(),
            Complex: HostGuestComplexKeyMaker(),
            NRotaxane: NRotaxaneKeyMaker(),
            Macrocycle: MacrocycleKeyMaker(),
        }

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return (
        )
