class FunctionalGroupFactory:
    """
    An abstract base class for functional group factories.

    The purpose of a functional group factory is to create
    :class:`.FunctionalGroup` instances. It allows the user to avoid
    creating :class:`.FunctionalGroup` instances manually. The user
    will use :meth:`.get_functional_groups` to achieve this. Subclasses
    of this class are made to customize the automatic creation
    process.

    Examples
    --------
    *Using a Single Functional Group From a Factory*

    You have a building block with two functional groups, but you want
    to use just one.

    .. testcode:: using-a-single-functional-group-factory

        import stk

        amino_factory = stk.PrimaryAminoFactory()
        building_block = stk.BuildingBlock('NCCN')
        amino_group1, *rest = amino_factory.get_functional_groups(
            molecule=building_block,
        )
        building_block = building_block.with_functional_groups(
            functional_groups=(amino_group1, ),
        )

    .. testcode:: using-a-single-functional-group-factory
        :hide:

        assert all(
            isinstance(functional_group, stk.PrimaryAmino)
            for functional_group
            in building_block.get_functional_groups()
        )
        assert building_block.get_num_functional_groups() == 1

    *Using a Subset of Functional Groups From a Factory*

    You have multiple functional groups, but you want the building
    block to use a specific subset.

    .. testcode:: using-a-subset-of-functional-groups-from-a-factory

        import stk

        bromo_factory = stk.BromoFactory()
        building_block = stk.BuildingBlock('BrCC(Br)CC(Br)CC(Br)CCBr')
        bromo_groups = tuple(bromo_factory.get_functional_groups(
            molecule=building_block,
        ))
        building_block = building_block.with_functional_groups(
            functional_groups=(
                bromo_groups[0],
                bromo_groups[3],
                bromo_groups[4],
            ),
        )

    .. testcode:: using-a-subset-of-functional-groups-from-a-factory
        :hide:

        assert all(
            isinstance(functional_group, stk.Bromo)
            for functional_group
            in building_block.get_functional_groups()
        )
        assert building_block.get_num_functional_groups() == 3

    More usage examples can be found in the docstrings of the
    various subclasses.

    *Subclass Implementation*

    The source of the subclasses, listed in
    :mod:`.functional_group_factory`, can serve as good examples.

    """

    def get_functional_groups(self, molecule):
        """
        Yield functional groups in `molecule`.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            The molecule, whose functional groups are to be found.

        Yields
        ------
        :class:`.FunctionalGroup`
            A functional group in `molecule`.

        Examples
        --------
        See :class:`.FunctionalGroupFactory`.

        """

        raise NotImplementedError()

    def __repr__(self) -> str:
        return f"{type(self).__name__}()"
