"""
Constructed Molecule Key Maker
==============================

"""

from ..molecule import InchiKey
from ..building_block import BuildingBlockKeyMaker


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
                building_blocks = tuple(
                    constructed_molecule.get_building_blocks()
                )
                return len(building_blocks)


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

    """

    def __init__(
        self,
        name='ConstructedMoleculeKey',
        molecule_key_maker=InchiKey(),
        building_block_key_maker=BuildingBlockKeyMaker(),
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

        """

        self._name = name
        self._molecule_key_maker = molecule_key_maker
        self._building_block_key_maker = building_block_key_maker
        self._topology_graph_key_makers = (
            self._get_default_topology_graph_key_makers()
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

    def __str__(self):
        return repr(self)

    def __repr__(self):
        name = (
            ''
            if self._name == 'ConstructedMoleculeKey'
            else f'name={self._name!r}'
        )
        molecule_key_maker = (
            ''
            if isinstance(self._molecule_key_maker, InchiKey)
            else f'molecule_key_maker={self._molecule_key_maker!r}'
        )
        building_block_key_maker = (
            ''
            if isinstance(
                self._building_block_key_maker,
                BuildingBlockKeyMaker,
            )
            else (
                f'building_block_key_maker='
                f'{self._building_block_key_maker!r}'
            )
        )
        parameters = (
            name,
            molecule_key_maker,
            building_block_key_maker,
        )
        parameters = ', '.join(
            parameter for parameter in parameters if parameter
        )
        return f'{self.__class__.__name__}({parameters})'
