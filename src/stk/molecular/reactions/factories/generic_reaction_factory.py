"""
Generic Reaction Factory
========================

"""

from ...functional_groups import (
    Aldehyde,
    Alkene,
    Alkyne,
    Amide,
    PrimaryAmino,
)
from ..reactions import OneOneReaction, OneTwoReaction, TwoTwoReaction
from .reaction_factory import ReactionFactory

# Impose the same interface on all reaction initializers.


def _one_one_reaction(
    construction_state,
    functional_group1,
    functional_group2,
    bond_order,
    periodicity
):
    return OneOneReaction(
        functional_group1=functional_group1,
        functional_group2=functional_group2,
        bond_order=bond_order,
        periodicity=periodicity,
    )


def _one_two_reaction(
    construction_state,
    functional_group1,
    functional_group2,
    bond_order,
    periodicity
):
    return OneTwoReaction(
        functional_group1=functional_group1,
        functional_group2=functional_group2,
        bond_order=bond_order,
        periodicity=periodicity,
    )


def _two_two_reaction(
    construction_state,
    functional_group1,
    functional_group2,
    bond_order,
    periodicity
):
    return TwoTwoReaction(
        construction_state=construction_state,
        functional_group1=functional_group1,
        functional_group2=functional_group2,
        bond_order=bond_order,
        periodicity=periodicity,
    )


def _get_reaction_key(functional_groups):
    """
    Return a key for :data:`._reactions`

    Parameters
    ----------
    functional_groups : :class:`iterable`
        An :class:`iterable` of :class:`.GenericFunctionalGroup`.
        The correct reaction must be selected for these functional
        groups.

    Returns
    -------
    :class:`.frozenset`
        A key for :data:`_reactions`, which maps to the correct
        reaction.

    """

    return frozenset(
        functional_group.get_num_bonders()
        for functional_group in functional_groups
    )


_reactions = {
    frozenset({1}): _one_one_reaction,
    frozenset({1, 2}): _one_two_reaction,
    frozenset({2}): _two_two_reaction,
}


class GenericReactionFactory(ReactionFactory):
    """
    Create reactions for :class:`.GenericFunctionalGroup` instances.

    This reaction factory assumes that the functional groups, which
    belong to the :class:`.EdgeGroup` passed to it, are
    :class:`.GenericFunctionalGroup` instances. It returns a
    :class:`.Reaction` suitable for two such instances.

    """

    def __init__(self, bond_orders=None):
        """
        Initialize a :class:`.GenericReactionFactory`.

        Parameters
        ----------
        bond_orders : :class:`dict`, optional
            Maps a :class:`frozenset` of
            :class:`.GenericFunctionalGroup` subclasses to the
            bond orders for their respective reactions, if a pair
            of functional groups is missing, a default bond order of
            1 will be used for their reactions. If `bond_orders` is
            ``None``, the following :class:`dict` will be used

            .. testcode:: init

                import stk

                bond_orders = {
                    frozenset({stk.PrimaryAmino, stk.Aldehyde}): 2,
                    frozenset({stk.PrimaryAmino, stk.Aldehyde}): 2,
                    frozenset({stk.Amide, stk.PrimaryAmino}): 2,
                    frozenset({stk.Alkene}): 2,
                    frozenset({stk.Alkyne}): 2,
                }

            This means that if you want to get a reaction for
            an amine and an aldehyde functional group, the reaction
            will create bonds with a bond order of 2.

        """

        # Used for __repr__.
        self._default_bond_orders = bond_orders is None

        if bond_orders is None:
            bond_orders = {
                frozenset({PrimaryAmino, Aldehyde}): 2,
                frozenset({Amide, Aldehyde}): 2,
                frozenset({Amide, PrimaryAmino}): 2,
                frozenset({Alkene}): 2,
                frozenset({Alkyne}): 2,
            }

        self._bond_orders = bond_orders

    def get_reaction(self, construction_state, edge_group):
        functional_groups = tuple(
            construction_state.get_edge_group_functional_groups(
                edge_group=edge_group,
            )
        )
        functional_group1, functional_group2 = functional_groups
        edge = construction_state.get_edge(
            edge_id=next(edge_group.get_edge_ids()),
        )
        return _reactions[_get_reaction_key(functional_groups)](
            construction_state=construction_state,
            functional_group1=functional_group1,
            functional_group2=functional_group2,
            bond_order=self._bond_orders.get(
                self._get_bond_order_key(functional_groups),
                1,
            ),
            periodicity=edge.get_periodicity(),
        )

    def _get_bond_order_key(self, functional_groups):
        """
        Return a key for the :attr:`_bond_orders`.

        Parameters
        ----------
        functional_groups : :class:`iterable`
            A :class:`iterable` of :class:`.GenericFunctionalGroup`
            instances. You want to to get the bond order for a reaction
            involving these instances.

        Returns
        -------
        :class:`frozenset`
            A key for :attr:`_bond_orders`, which maps to the correct
            bond order.

        """

        return frozenset(map(type, functional_groups))

    def __str__(self):
        return repr(self)

    def __repr__(self):
        bond_orders = (
            '' if self._default_bond_orders else f'{self._bond_orders}'
        )
        return f'{self.__class__.__name__}({bond_orders})'
