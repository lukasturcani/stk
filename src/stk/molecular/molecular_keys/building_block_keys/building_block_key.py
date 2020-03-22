"""
Building Block Key
==================

"""

from ..molecule_keys import InchiKey
from .functional_group_keys import (
    AlcoholKey,
    AldehydeKey,
    AlkeneKey,
    AlkyneKey,
    AmideKey,
    BoronicAcidKey,
    BromoKey,
    CarboxylicAcidKey,
    DibromoKey,
    DifluoroKey,
    DiolKey,
    FluoroKey,
    FunctionalGroupKey,
    GenericFunctionalGroupKey,
    IodoKey,
    PrimaryAminoKey,
    RingAmineKey,
    SecondaryAminoKey,
    ThioacidKey,
    ThiolKey,
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


class BuildingBlockKey:
    """
    An abstract base class for :class:`.BuildingBlock` keys.

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

    There's just two methods two simple methods to implement.

    .. code-block:: python

        import stk

        class NumFunctionalGroups(stk.BuildingBlockKey):
            def name(self):
                # What string this is completely up to you. It does
                # not have to be related to the class name.
                return 'num_functional_groups'

            def get_key(self, building_block):
                return building_block.get_num_functional_groups()

        # A usage example of the new subclass.

        jsonizer = stk.BuildingBlockJsonizer(
            building_block_keys=(NumFunctionalGroups(), ),
        )
        # Create a JSON representation of a building block, which
        # holds the number of functional groups.
        json = jsonizer.to_json(stk.BuildingBlock('NCCN'))

    *Usage*

    Because :class:`.BuildingBlockKey` comes with a default
    implementation, it can be used directly, instead of having to
    make a subclass

    .. code-block:: python

        import stk

        jsonizer = stk.BuildingBlockJsonizer(
            building_block_keys=(stk.BuildingBlockKey(), ),
        )
        # Create a JSON representation of a building block, which
        # holds the default building block key.
        json = jsonizer.to_json(stk.BuildingBlock('NCCN'))

    You want to use an InChI instead of an InChIKey in the building
    block key

    .. code-block:: python

        import stk
        jsonizer = stk.BuildingBlockJsonizer(
            building_block_keys=(
                stk.BuildingBlockKey(
                    # If you change the nature of the key, its a good
                    # idea to change its name to reflect that.
                    name='InChIBuildingBlockKey',
                    molecule_key=stk.Inchi(),
                ),
                # You can still keep the default building block key
                # too. No pressure though, excluding it from
                # this tuple is also valid.
                stk.BuildingBlockKey(),
            ),
        )
        # Create a JSON representation of a building block, which
        # holds two building block keys, one featuring an InChI and
        # one featuring an InChIKey.
        json = jsonizer.to_json(stk.BuildingBlock('NCCN'))

    You want to use a different functional group key for some
    :class:`.FunctionalGroup` subclass

    .. code-block:: python

        import stk

        # Define the alternate functional group key you want to use.
        # In this case, the bromo key is given by the id of
        # the bromine atom.
        class MyBromoKey:
            def get_key(self, bromo):
                return f'BromoKey({bromo.get_bromine().get_id()}'


        building_block_key = stk.BuildingBlockKey(
            # Use the alternate key instead of the default.
            functional_group_keys={
                stk.Bromo: MyBromoKey(),
            },
        )

        # Use building_block_key as normal.
        jsonizer = stk.BuildingBlockJsonizer(
            building_block_keys=(building_block_key, ),
        )
        json = jsonizer.to_json(
            building_block=stk.BuildingBlock(
                smiles='BrCC',
                functional_groups=[stk.BromoFactory()],
            ),
        )

    If you want to use a key for a new :class:`.FunctionalGroup` the
    process is essentially the same the one described above

    .. code-block:: python

        # Define your functional group subclass
        class MyFunctionalGroup(stk.FunctionalGroup):
            # Your implementation goes here.


        # Define the functional group key for your subclass
        class MyFunctionalGroupKey:
            def get_key(self, my_functional_group):
                # Get the key of my_functional_group somehow.


        building_block_key = stk.BuildingBlockKey(
            # Use the updated mapping.
            functional_group_keys={
                MyFunctionalGroup: MyFunctionalGroupKey(),
            },
        )

    """

    def __init__(
        self,
        name='BuildingBlockKey',
        molecule_key=InchiKey(),
        functional_group_keys=None,
    ):
        """
        Initialize a :class:`.BuildingBlockKey` instance.

        Parameters
        ----------
        name : :class:`str`, optional
            The name of the key.

        molecule_key : :class:`.MoleculeKey`, optional
            Used to generate the part of key responsible for the
            molecular component of a :class:`.BuildingBlock`.

        functional_group_keys : :class:`dict`, optional
            Map a :class:`FunctionalGroup` subclass to the
            functional group key, which should be used to get the
            key for its instances. For a list of built-in
            functional group keys, see :mod:`.functional_group_keys`.
            If ``None``, the built-in keys will be used.

            This parameter will only update the default mapping.
            Providing an empty ``{}``, will mean that all the defaults
            are kept. Providing

            .. code-block:: python

                building_block_key = stk.BuildingBlockKey(
                    functional_group_keys={
                        stk.Bromo: CustomBromoKey(),
                    }
                )

            will keep all the defaults, with the exception of
            :class:`.Bromo`, which will use :class:`.CustomBromoKey`
            instead of the default.

        """

        if functional_group_keys is None:
            functional_group_keys = {}

        self._name = name
        self._molecule_key = molecule_key
        self._functional_group_keys = (
            self._get_default_functional_group_keys()
        )
        self._functional_group_keys.update(functional_group_keys)

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

        functional_group_keys = '-'.join(
            str(self._functional_group_keys[type(fg)].get_key(fg))
            for fg in building_block.get_functional_groups()
        )
        placer_ids = ''.join(
            str(id_) for id_ in building_block.get_placer_ids()
        )
        return (
            f'{self._molecule_key.get_key(building_block)}-'
            f'{functional_group_keys}-'
            f'{placer_ids}'
        )

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._name}, {self._molecule_key!r}, '
            f'{self._functional_group_keys!r}'
            ')'
        )

    @staticmethod
    def _get_default_functional_group_keys():
        """
        Get the default functional group key for each functional group.

        Returns
        -------
        :class:`dict`
            Maps a :class:`.FunctionalGroup` subclass to the
            :mod:`functional_group_key <.functional_group_keys>` that
            should be used with it.

        """

        return {
            Alcohol: AlcoholKey(),
            Aldehyde: AldehydeKey(),
            Alkene: AlkeneKey(),
            Alkyne: AlkyneKey(),
            Amide: AmideKey(),
            BoronicAcid: BoronicAcidKey(),
            Bromo: BromoKey(),
            CarboxylicAcid: CarboxylicAcidKey(),
            Dibromo: DibromoKey(),
            Difluoro: DifluoroKey(),
            Diol: DiolKey(),
            Fluoro: FluoroKey(),
            FunctionalGroup: FunctionalGroupKey(),
            GenericFunctionalGroup: GenericFunctionalGroupKey(),
            Iodo: IodoKey(),
            PrimaryAmino: PrimaryAminoKey(),
            RingAmine: RingAmineKey(),
            SecondaryAmino: SecondaryAminoKey(),
            Thioacid: ThioacidKey(),
            Thiol: ThiolKey(),
        }
