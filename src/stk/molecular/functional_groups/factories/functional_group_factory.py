"""
Functional Group Factory
========================

#. :class:`.AlcoholFactory`
#. :class:`.AldehydeFactory`
#. :class:`.AmideFactory`
#. :class:`.BoronicAcidFactory`
#. :class:`.BromoFactory`
#. :class:`.CarboxylicAcidFactory`
#. :class:`.DibromoFactory`
#. :class:`.DifluoroFactory`
#. :class:`.DiolFactory`
#. :class:`.FluoroFactory`
#. :class:`.IodoFactory`
#. :class:`.PrimaryAminoFactory`
#. :class:`.RingAmineFactory`
#. :class:`.SecondaryAminoFactory`
#. :class:`.TerminalAlkeneFactory`
#. :class:`.TerminalAlkyneFactory`
#. :class:`.ThioacidFactory`
#. :class:`.ThiolFactory`

"""


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
    *Usage*

    You have a building block with two functional groups, but you want
    to use just one.

    .. code-block:: python

        import stk

        amino_factory = stk.PrimaryAminoFactory()
        building_block = stk.BuildingBlock('NCCN')
        amino_group = next(
            amino_factory.get_functional_groups(building_block)
        )
        building_block = building_block.with_functional_groups(
            functional_groups=(amino_group, ),
        )

    You have multiple functional groups, but you want the building
    block to use a specific subset.

    .. code-block:: python

        import stk

        bromo_factory = stk.BromoFactory()
        building_block = stk.BuildingBlock('BrCC(Br)CC(Br)CC(Br)CCBr')
        bromo_groups = tuple(
            bromo_factory.get_functional_groups(molecule)
        )
        building_block = building_block.with_functional_groups(
            functional_groups=(
                bromo_groups[0],
                bromo_groups[3],
                bromo_groups[4],
            ),
        )

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
