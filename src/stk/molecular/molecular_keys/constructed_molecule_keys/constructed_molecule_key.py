"""
Constructed Molecule Key
========================

"""

from .topology_graph_keys import (
    CageKey,
    CofKey,
    LinearPolymerKey,
    HostGuestComplexKey,
    NRotaxaneKey,
    MacrocycleKey,
)
from ..molecule_keys import InchiKey
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


class ConstructedMoleculeKey:
    """
    An abstract base class for :class:`.ConstructedMolecule` keys.

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


        class NumBuildingBlocks(stk.ConstructedMoleculeKey):
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
            constructed_molecule_keys=(NumBuildingBlocks(), ),
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

    Because :class:`.ConstructedMoleculeKey` comes with a default
    implementation, it can be used directly, instead of having to
    make a subclass

    .. code-block:: python

        import stk

        jsonizer = stk.ConstructedMoleculeJsonizer(
            constructed_molecule_keys=(stk.ConstructedMoleculeKey(), ),
        )
        # Create a JSON representation of a constructed molecule, which
        # holds the default constructed molecule key.
        json = jsonizer.to_json(polymer)

    You want to use an InChI instead of an InChIKey in the
    constructed molecule key

    .. code-block:: python

        jsonizer = stk.ConstructedMoleculeJsonizer(
            constructed_molecule_keys=(
                stk.ConstructedMoleculeKey(
                    # If you change the nature of the key, its a good
                    # idea to change its name to reflect that.
                    name='InChIConstructedMoleculeKey',
                    molecule_key=stk.Inchi(),
                ),
                # You can still keep the default constructed molecule
                # key too. No pressure though, excluding it from
                # this tuple is also valid.
                stk.ConstructedMoleculeKey(),
            ),
        )
        # Create a JSON representation of a constructed molecule, which
        # holds two constructed molecule keys, one featuring an
        # InChI and one featuring an InChIKey.
        json = jsonizer.to_json(polymer)

    You want to use a different topology graph key for some
    :class:`.TopologyGraph` subclass

    .. code-block:: python

        import stk

        # Define the alternate topology graph key you want to use.
        # In this case, the polymer key is given by the number of
        # building blocks.
        class LinearKey:
            def get_key(self, linear):
                building_blocks = (
                    linear.get_building_blocks()
                )
                return sum(1 for _ in building_blocks)

        constructed_molecule_key = stk.ConstructedMoleculeKey(
            # Use an alternate topology graph key for
            # stk.polymer.Linear instances.
            topology_graph_keys={
                stk.polymer.Linear: LinearKey(),
            }
        )

        # Use the constructed_molecule_key as normal.
        jsonizer = stk.ConstructedMoleculeJsonizer(
            constructed_molecule_keys=(
                constructed_molecule_key,
            ),
        )
        json = jsonizer.to_json(polymer)

    If you want to use a key for a new :class:`.TopologyGraph` the
    process is essentially the same as the one described above

    .. code-block:: python

        # Define your TopologyGraph subclass
        class MyTopologyGraph(stk.TopologyGraph):
            # Your implementation goes here


        # Define the topology graph key for your subclass
        class MyTopologyGraphKey:
            def get_key(self, my_topology_graph):
                # Get the key for my_topology_graph somehow.


        constructed_molecule_key = stk.ConstructedMoleculeKey(
            # use the new key.
            topology_graph_keys={
                MyTopologyGraph: MyTopologyGraphKey(),
            },
        )

    """

    def __init__(
        self,
        name='ConstructedMoleculeKey',
        molecule_key=InchiKey(),
        topology_graph_keys=None,
    ):
        """
        Initialize a :class:`.ConstructedMoleculeKey` instance.

        Parameters
        ----------
        name : :class:`str`, optional
            The name of the key.

        molecule_key : :class:`.MoleculeKey`, optional
            Used to generate the part of the key responsible for the
            molecular component of :class:`.ConstructedMolecule`.

        topology_graph_keys : :class:`dict`, optional
            Maps a :class:`.TopologyGraph` subclass to the
            topology graph key, which should be used to get key for
            its instances. For a list of built-in topology graph
            keys, see :mod:`.topology_graph_keys`. If ``None``,
            the built-in keys will be used.

            This parameter will only update the default mapping.
            Providing an empty ``{}``, will mean that all the
            defaults are kept. Providing

            .. code-block:: python

                constructed_molecule_key = stk.ConstructedMOleculeKey(
                    topology_graph_keys={
                        stk.polymer.Linear: CustomLinearKey(),
                    },
                )

            will keep all the defaults, with the exception of
            :class:`.polymer.Linear`, which will use
            :class:`.CustomLinearKey` instead of the default.

        """

        if topology_graph_keys is None:
            topology_graph_keys = {}

        self._name = name
        self._molecule_key = molecule_key
        self._topology_graph_keys = (
            self._get_default_topology_graph_keys()
        )
        self._topology_graph_key.update(topology_graph_keys)
        # Used for __repr__().
        self._input_topology_graph_keys = topology_graph_keys

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
        key = type(self._topology_graph_keys[topology_graph])
        return key.get_key(topology_graph)

    @staticmethod
    def _get_default_topology_graph_keys():
        """
        Get the default topology graph key for each topology graph.

        Returns
        -------
        :class:`dict`

        """

        return {
            EightPlusSixteen: CageKey(),
            EightPlusTwelve: CageKey(),
            FivePlusTen: CageKey(),
            FourPlusEight: CageKey(),
            FourPlusFour: CageKey(),
            FourPlusSix: CageKey(),
            FourPlusSix2: CageKey(),
            OnePlusOne: CageKey(),
            SixPlusEight: CageKey(),
            SixPlusNine: CageKey(),
            SixPlusTwelve: CageKey(),
            TenPlusTwenty: CageKey(),
            ThreePlusSix: CageKey(),
            TwelvePlusThirty: CageKey(),
            TwentyPlusThirty: CageKey(),
            TwoPlusFour: CageKey(),
            TwoPlusThree: CageKey(),
            TwoPlusTwo: CageKey(),
            Hexagonal: CofKey(),
            Honeycomb: CofKey(),
            Kagome: CofKey(),
            LinkerlessHoneycomb: CofKey(),
            Square: CofKey(),
            Linear: LinearPolymerKey(),
            Complex: HostGuestComplexKey(),
            NRotaxane: NRotaxaneKey(),
            Macrocycle: MacrocycleKey(),
        }

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return (
        )
