"""
Building Block Key
==================

"""


class BuildingBlockKey:
    """
    An abstract base class for building block keys.

    Notes
    -----
    You might notice that the public methods of this abstract base
    class are implemented. This is purely for convenience when
    implementing subclasses. The implemented public methods are
    simply default implementations, which can be safely ignored or
    overridden, when implementing subclasses.

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

        num_functional_groups = stk.BuildingBlockKey(
            name='num_functional_groups',
            key=lambda building_block:
                building_block.get_num_functional_groups(),
        )
        jsonizer = stk.BuildingBlockJsonizer(
            building_block_keys=(num_functional_groups, ),
        )
        # Create a JSON representation of a building block, which
        # holds the number of functional groups.
        json = jsonizer.to_json(stk.BuildingBlock('NCCN'))

    """

    def __init__(self, molecule_key, functional_group_key):
        """

        """

        self._molecule_key = molecule_key
        self._functional_group_key = functional_group_key

    def get_name(self):
        """
        """

        return self._name

    def get_key(self, building_block):
        """

        """

        functional_group_keys = '-'.join(
            self._functional_group_key(functional_group)
            for functional_group
            in building_block.get_functional_groups()
        )
        placer_ids = ''.join(
            str(id_) for id_ in building_block.get_placer_ids()
        )
        return (
            f'{self._molecule_key(building_block)}-'
            f'{functional_group_keys}-'
            f'{placer_ids}'
        )
