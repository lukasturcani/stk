"""
Building Block Key Maker
========================

"""

from ..molecule import InchiKey
from .functional_group import (
    AlcoholKeyMaker,
    AldehydeKeyMaker,
    AlkeneKeyMaker,
    AlkyneKeyMaker,
    AmideKeyMaker,
    BoronicAcidKeyMaker,
    BromoKeyMaker,
    CarboxylicAcidKeyMaker,
    DibromoKeyMaker,
    DifluoroKeyMaker,
    DiolKeyMaker,
    FluoroKeyMaker,
    FunctionalGroupKeyMaker,
    GenericFunctionalGroupKeyMaker,
    IodoKeyMaker,
    PrimaryAminoKeyMaker,
    RingAmineKeyMaker,
    SecondaryAminoKeyMaker,
    ThioacidKeyMaker,
    ThiolKeyMaker,
)
from ...functional_groups import (
    Alcohol,
    Aldehyde,
    Alkene,
    Alkyne,
    Amide,
    BoronicAcid,
    Bromo,
    CarboxylicAcid,
    Dibromo,
    Difluoro,
    Diol,
    Fluoro,
    FunctionalGroup,
    GenericFunctionalGroup,
    Iodo,
    PrimaryAmino,
    RingAmine,
    SecondaryAmino,
    Thioacid,
    Thiol,
)


class BuildingBlockKeyMaker:
    """
    An abstract base class for :class:`.BuildingBlock` key makers.

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


        class NumFunctionalGroups(stk.BuildingBlockKeyMaker):
            def __init__(self):
                return

            def name(self):
                # What string this is is completely up to you. It does
                # not have to be related to the class name.
                return 'num_functional_groups'

            def get_key(self, building_block):
                return building_block.get_num_functional_groups()


        # A usage example of the new subclass.

        jsonizer = stk.BuildingBlockJsonizer(
            building_block_key_makers=(NumFunctionalGroups(), ),
        )
        # Create a JSON representation of a building block, which
        # holds the number of functional groups.
        json = jsonizer.to_json(stk.BuildingBlock('NCCN'))

    *Usage*

    Because :class:`.BuildingBlockKeyMaker` comes with a default
    implementation, it can be used directly, instead of having to
    make a subclass

    .. code-block:: python

        import stk

        jsonizer = stk.BuildingBlockJsonizer(
            building_block_key_makers=(stk.BuildingBlockKeyMaker(), ),
        )
        # Create a JSON representation of a building block, which
        # holds the default building block key.
        json = jsonizer.to_json(stk.BuildingBlock('NCCN'))

    You want to use an InChI instead of an InChIKey in the building
    block key maker

    .. code-block:: python

        import stk
        jsonizer = stk.BuildingBlockJsonizer(
            building_block_key_makers=(
                stk.BuildingBlockKeyMaker(
                    # If you change the nature of the key, its a good
                    # idea to change its name to reflect that.
                    name='InChIBuildingBlockKey',
                    molecule_key_maker=stk.Inchi(),
                ),
                # You can still keep the default building block key
                # maker too. No pressure though, excluding it from
                # this tuple is also valid.
                stk.BuildingBlockKeyMaker(),
            ),
        )
        # Create a JSON representation of a building block, which
        # holds two building block keys, one featuring an InChI and
        # one featuring an InChIKey.
        json = jsonizer.to_json(stk.BuildingBlock('NCCN'))

    You want to use a different functional group key maker for some
    :class:`.FunctionalGroup` subclass

    .. code-block:: python

        import stk

        # Define the alternate functional group key maker you want to
        # use. In this case, the bromo key is given by the id of
        # the bromine atom.
        class MyBromoKeyMaker:
            def get_key(self, bromo):
                return f'BromoKey({bromo.get_bromine().get_id()}'


        building_block_key_maker = stk.BuildingBlockKeyMaker(
            # Use the alternate key maker instead of the default.
            functional_group_key_makers={
                stk.Bromo: MyBromoKeyMaker(),
            },
        )

        # Use building_block_key_maker as normal.
        jsonizer = stk.BuildingBlockJsonizer(
            building_block_key_makers=(building_block_key_maker, ),
        )
        json = jsonizer.to_json(
            building_block=stk.BuildingBlock(
                smiles='BrCC',
                functional_groups=[stk.BromoFactory()],
            ),
        )

    If you want to use a key maker for a new :class:`.FunctionalGroup`
    the process is essentially the same the one described above

    .. code-block:: python

        # Define your functional group subclass
        class MyFunctionalGroup(stk.FunctionalGroup):
            # Your implementation goes here.


        # Define the functional group key maker for your subclass
        class MyFunctionalGroupKeyMaker:
            def get_key(self, my_functional_group):
                # Get the key of my_functional_group somehow.


        building_block_key_maker = stk.BuildingBlockKeyMaker(
            # Use the new key maker.
            functional_group_key_makers={
                MyFunctionalGroup: MyFunctionalGroupKeyMaker(),
            },
        )

    """

    def __init__(
        self,
        name='BuildingBlockKey',
        molecule_key_maker=InchiKey(),
        functional_group_key_makers=None,
    ):
        """
        Initialize a :class:`.BuildingBlockKeyMaker` instance.

        Parameters
        ----------
        name : :class:`str`, optional
            The name of the key made by the maker.

        molecule_key_maker : :class:`.MoleculeKeyMaker`, optional
            Used to generate the part of key responsible for the
            molecular component of a :class:`.BuildingBlock`.

        functional_group_key_makers : :class:`dict`, optional
            Maps a :class:`.FunctionalGroup` subclass to the
            functional group key maker, which should be used to get the
            key for its instances. For a list of built-in
            functional group key makers, see
            :mod:`.functional_group_key_makers`.
            If ``None``, the built-in key makers will be used.

            This parameter will only update the default mapping.
            Providing an empty ``{}``, will mean that all the defaults
            are kept. Providing

            .. code-block:: python

                building_block_key_maker = stk.BuildingBlockKeyMaker(
                    functional_group_key_makers={
                        stk.Bromo: CustomBromoKeyMaker(),
                    }
                )

            will keep all the defaults, with the exception of
            :class:`.Bromo`, which will use
            :class:`.CustomBromoKeyMaker` instead of the default.

        """

        if functional_group_key_makers is None:
            functional_group_key_makers = {}

        self._name = name
        self._molecule_key_maker = molecule_key_maker
        self._functional_group_key_makers = (
            self._get_default_functional_group_key_makers()
        )
        self._functional_group_key_makers.update(
            functional_group_key_makers
        )
        # Used for __repr__().
        self._input_functional_group_key_makers = (
            functional_group_key_makers
        )

    def get_name(self):
        """
        Get the name of the key.

        Returns
        -------
        :class:`str`
            The name of the key.

        """

        return self._name

    def get_key(self, building_block):
        """
        Get the key of `building_block`.

        Parameters
        ----------
        building_block : :class:`.BuildingBlock`
            The building block for which a key is wanted.

        Returns
        -------
        :class:`object`
            The key of `building_block`.

        """

        building_block = building_block.with_canonical_atom_ordering()

        def atom_ids(functional_group):
            return sorted(functional_group.get_atom_ids())

        functional_groups = sorted(
            building_block.get_functional_groups(),
            key=atom_ids,
        )
        key_makers = map(
            self._functional_group_key_makers.get,
            map(type, functional_groups),
        )
        functional_group_keys = '-'.join(
            str(key_maker.get_key(fg))
            for key_maker, fg in zip(key_makers, functional_groups)
        )
        placer_ids = ''.join(
            str(id_) for id_ in sorted(building_block.get_placer_ids())
        )
        return (
            f'{self._molecule_key_maker.get_key(building_block)}-'
            f'{functional_group_keys}-'
            f'{placer_ids}'
        )

    def __str__(self):
        return repr(self)

    def __repr__(self):
        name = (
            ''
            if self._name == 'BuildingBlockKey'
            else f'name={self._name!r}'
        )
        molecule_key_maker = (
            ''
            if isinstance(self._molecule_key_maker, InchiKey)
            else f'molecule_key_maker={self._molecule_key_maker!r}'
        )
        functional_group_key_makers = (
            ''
            if not self._input_functional_group_key_makers
            else (
                'functional_group_key_makers='
                f'{self._input_functional_group_key_makers!r}'
            )
        )
        return (
            f'{self.__class__.__name__}('
            f'{name}'
            f'{", " if molecule_key_maker else ""}'
            f'{molecule_key_maker}'
            f'{", " if functional_group_key_makers else ""}'
            f'{functional_group_key_makers}'
            ')'
        )

    @staticmethod
    def _get_default_functional_group_key_makers():
        """
        Get the default functional group key makers.

        Returns
        -------
        :class:`dict`
            Maps a :class:`.FunctionalGroup` subclass to the
            functional group key maker that should be used with it.

        """

        return {
            Alcohol: AlcoholKeyMaker(),
            Aldehyde: AldehydeKeyMaker(),
            Alkene: AlkeneKeyMaker(),
            Alkyne: AlkyneKeyMaker(),
            Amide: AmideKeyMaker(),
            BoronicAcid: BoronicAcidKeyMaker(),
            Bromo: BromoKeyMaker(),
            CarboxylicAcid: CarboxylicAcidKeyMaker(),
            Dibromo: DibromoKeyMaker(),
            Difluoro: DifluoroKeyMaker(),
            Diol: DiolKeyMaker(),
            Fluoro: FluoroKeyMaker(),
            FunctionalGroup: FunctionalGroupKeyMaker(),
            GenericFunctionalGroup: GenericFunctionalGroupKeyMaker(),
            Iodo: IodoKeyMaker(),
            PrimaryAmino: PrimaryAminoKeyMaker(),
            RingAmine: RingAmineKeyMaker(),
            SecondaryAmino: SecondaryAminoKeyMaker(),
            Thioacid: ThioacidKeyMaker(),
            Thiol: ThiolKeyMaker(),
        }
