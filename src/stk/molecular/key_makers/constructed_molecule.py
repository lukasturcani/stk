"""
Constructed Molecule Key Maker
==============================

#. :class:`.ConstructedMoleculeKeyMaker`
#. :class:`.Inchi`
#. :class:`.InchiKey`
#. :class:`.MoleculeKeyMaker`

"""

from .molecule import InchiKey


class ConstructedMoleculeKeyMaker:
    """
    Abstract base class for making :class:`.ConstructedMolecule` keys.

    Note that every :class:`.MoleculeKeyMaker` is also a valid
    :class:`.ConstructedMoleculeKeyMaker`.

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

            def key_name(self):
                # What string this is is completely up to you. It does
                # not have to be related to the class name.
                return 'num_building_blocks'

            def get_key(self, molecule):
                building_blocks = tuple(molecule.get_building_blocks())
                return len(building_blocks)


        # A usage example of the new subclass.

        jsonizer = stk.ConstructedMoleculeJsonizer(
            key_makers=(NumBuildingBlocks(), ),
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
            key_makers=(
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
            key_makers=(
                stk.ConstructedMoleculeKeyMaker(
                    # If you change the nature of the key, its a good
                    # idea to change its name to reflect that.
                    key_name='InChIConstructedMoleculeKey',
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
        key_name='ConstructedMoleculeKey',
        molecule_key_maker=InchiKey(),
    ):
        """
        Initialize a :class:`.ConstructedMoleculeKeyMaker` instance.

        Parameters
        ----------
        key_name : :class:`str`, optional
            The name of the key made by the maker.

        molecule_key_maker : :class:`.MoleculeKeyMaker`, optional
            Used to generate the part of the key responsible for the
            molecular component of :class:`.ConstructedMolecule`.

        """

        self._key_name = key_name
        self._molecule_key_maker = molecule_key_maker

    def get_key_name(self):
        """
        Get the name of the key.

        Returns
        -------
        :class:`str`
            The name.

        """

        return self._key_name

    def get_key(self, molecule):
        """
        Get the key of `molecule`.

        Parameters
        ----------
        molecule : :class:`.ConstructedMolecule`
            The constructed molecule for which a key is wanted.

        Returns
        -------
        :class:`object`

        """

        molecule_key = self._molecule_key_maker.get_key(molecule)
        building_blocks = tuple(map(
            self._molecule_key_maker.get_key,
            molecule.get_building_blocks(),
        ))
        num_building_blocks = tuple(map(
            molecule.get_num_building_block,
            molecule.get_building_blocks(),
        ))
        return str((
            molecule_key,
            building_blocks,
            num_building_blocks,
        ))

    def __str__(self):
        return repr(self)

    def __repr__(self):
        key_name = (
            ''
            if self._key_name == 'ConstructedMoleculeKey'
            else f'key_name={self._key_name!r}'
        )
        molecule_key_maker = (
            ''
            if isinstance(self._molecule_key_maker, InchiKey)
            else f'molecule_key_maker={self._molecule_key_maker!r}'
        )
        parameters = ', '.join(
            parameter
            for parameter in (key_name, molecule_key_maker)
            if parameter
        )
        return f'{self.__class__.__name__}({parameters})'
