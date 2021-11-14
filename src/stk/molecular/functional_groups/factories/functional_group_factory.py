"""
Functional Group Factory
========================

.. toctree::
    :maxdepth: 2

    Alcohol Factory <\
stk.molecular.functional_groups.factories.alcohol_factory\
>
    Aldehyde Factory <\
stk.molecular.functional_groups.factories.aldehyde_factory\
>
    Amide Factory <\
stk.molecular.functional_groups.factories.amide_factory\
>
    Boronic Acid Factory <\
stk.molecular.functional_groups.factories.boronic_acid_factory\
>
    Bromo Factory <\
stk.molecular.functional_groups.factories.bromo_factory\
>
    Carboxylic Acid Factory <\
stk.molecular.functional_groups.factories.carboxylic_acid_factory\
>
    Dibromo Factory <\
stk.molecular.functional_groups.factories.dibromo_factory\
>
    Difluoro Factory <\
stk.molecular.functional_groups.factories.difluoro_factory\
>
    Diol Factory <\
stk.molecular.functional_groups.factories.diol_factory\
>
    Fluoro Factory <\
stk.molecular.functional_groups.factories.fluoro_factory\
>
    Iodo Factory <\
stk.molecular.functional_groups.factories.iodo_factory\
>
    Primary Amino Factory <\
stk.molecular.functional_groups.factories.primary_amino_factory\
>
    Ring Amine Factory <\
stk.molecular.functional_groups.factories.ring_amine_factory\
>
    Secondary Amino Factory <\
stk.molecular.functional_groups.factories.secondary_amino_factory\
>
    SMARTS Functional Group Factory <\
stk.molecular.functional_groups.factories.smarts_functional_group_factory\
>
    Terminal Alkene Factory <\
stk.molecular.functional_groups.factories.terminal_alkene_factory\
>
    Terminal Alkyne Factory <\
stk.molecular.functional_groups.factories.terminal_alkyne_factory\
>
    Thioacid Factory <\
stk.molecular.functional_groups.factories.thioacid_factory\
>
    Thiol Factory <\
stk.molecular.functional_groups.factories.thiol_factory\
>

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
