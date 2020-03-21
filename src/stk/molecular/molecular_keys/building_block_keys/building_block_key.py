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
    GenericFunctionalGroupKey,
    IodoKey,
    PrimaryAminoKey,
    RingAmineKey,
    SecondaryAminoKey,
    ThioAcidKey,
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
    GenericFunctionalGroup,
    Iodo,
    PrimaryAmino,
    RingAmine,
    SecondaryAmino,
    ThioAcid,
    Thiol,
)


class BuildingBlockKey:
    """
    An abstract base class for building block keys.

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

    You want to use a different functional group key for some
    :class:`.FunctionalGroup` subclass.

    . code-block:: python

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

        """

        if functional_group_keys is None:
            functional_group_keys = (
                self.get_default_functional_group_keys()
            )

        self._name = name
        self._molecule_key = molecule_key
        self._functional_group_keys = functional_group_keys

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
            self._functional_group_key.get_key(functional_group)
            for functional_group
            in building_block.get_functional_groups()
        )
        placer_ids = ''.join(
            str(id_) for id_ in building_block.get_placer_ids()
        )
        return (
            f'{self._molecule_key.get_key(building_block)}-'
            f'{functional_group_keys}-'
            f'{placer_ids}'
        )

    @staticmethod
    def get_default_functional_group_keys():
        """

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
            GenericFunctionalGroup: GenericFunctionalGroupKey(),
            Iodo: IodoKey(),
            PrimaryAmino: PrimaryAminoKey(),
            RingAmine: RingAmineKey(),
            SecondaryAmino: SecondaryAminoKey(),
            ThioAcid: ThioAcidKey(),
            Thiol: ThiolKey(),
        }
