"""
Generic Reaction Factory
========================

"""


import typing
from collections import abc

from .reaction_factory import ReactionFactory
from ..reactions import (
    Reaction,
    OneOneReaction,
    OneTwoReaction,
    TwoTwoReaction,
)
from ...topology_graphs import ConstructionState, EdgeGroup
from ...functional_groups import (
    GenericFunctionalGroup,
    Alkene,
    Alkyne,
    PrimaryAmino,
    Aldehyde,
    Amide,
)

__all__ = (
    'GenericReactionFactory',
)

# Impose the same interface on all reaction initializers.


def _one_one_reaction(
    construction_state: ConstructionState,
    functional_group1: GenericFunctionalGroup,
    functional_group2: GenericFunctionalGroup,
    bond_order: int,
    periodicity: tuple[int, int, int],
) -> OneOneReaction:
    return OneOneReaction(
        functional_group1=functional_group1,
        functional_group2=functional_group2,
        bond_order=bond_order,
        periodicity=periodicity,
    )


def _one_two_reaction(
    construction_state: ConstructionState,
    functional_group1: GenericFunctionalGroup,
    functional_group2: GenericFunctionalGroup,
    bond_order: int,
    periodicity: tuple[int, int, int],
) -> OneTwoReaction:
    return OneTwoReaction(
        functional_group1=functional_group1,
        functional_group2=functional_group2,
        bond_order=bond_order,
        periodicity=periodicity,
    )


def _two_two_reaction(
    construction_state: ConstructionState,
    functional_group1: GenericFunctionalGroup,
    functional_group2: GenericFunctionalGroup,
    bond_order: int,
    periodicity: tuple[int, int, int],
) -> TwoTwoReaction:
    return TwoTwoReaction(
        construction_state=construction_state,
        functional_group1=functional_group1,
        functional_group2=functional_group2,
        bond_order=bond_order,
        periodicity=periodicity,
    )


def _get_reaction_key(
    functional_groups: abc.Iterable[GenericFunctionalGroup],
) -> frozenset[int]:
    """
    Return a key for :data:`._reactions`

    Parameters:

        functional_groups:
            An :class:`iterable` of :class:`.GenericFunctionalGroup`.
            The correct reaction must be selected for these functional
            groups.

    Returns:

        A key for :data:`_reactions`, which maps to the correct
        reaction.

    """

    return frozenset(
        functional_group.get_num_bonders()
        for functional_group in functional_groups
    )


_GenericFunctionalGroups = frozenset[
    typing.Type[GenericFunctionalGroup]
]


class _GetReaction(typing.Protocol):
    def __call__(
        self,
        construction_state: ConstructionState,
        functional_group1: GenericFunctionalGroup,
        functional_group2: GenericFunctionalGroup,
        bond_order: int,
        periodicity: tuple[int, int, int],
    ) -> Reaction:

        pass


_reactions: dict[frozenset[int], _GetReaction] = {
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

    def __init__(
        self,
        bond_orders:
            typing.Optional[dict[_GenericFunctionalGroups, int]]
            = None,
    ) -> None:
        """
        Initialize a :class:`.GenericReactionFactory`.

        Parameters:

            bond_orders:
                Maps a :class:`frozenset` of
                :class:`.GenericFunctionalGroup` subclasses to the
                bond orders for their respective reactions, if a pair
                of functional groups is missing, a default bond order
                of 1 will be used for their reactions. If `bond_orders`
                is ``None``, the following :class:`dict` will be used

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

    def get_reaction(
        self,
        construction_state: ConstructionState,
        edge_group: EdgeGroup,
    ) -> Reaction:

        functional_groups = typing.cast(
            tuple[GenericFunctionalGroup, ...],
            tuple(
                construction_state.get_edge_group_functional_groups(
                    edge_group=edge_group,
                )
            ),
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

    def _get_bond_order_key(
        self,
        functional_groups: abc.Iterable[GenericFunctionalGroup],
    ) -> _GenericFunctionalGroups:
        """
        Return a key for the :attr:`_bond_orders`.

        Parameters:

            functional_groups:
                A :class:`iterable` of :class:`.GenericFunctionalGroup`
                instances. You want to to get the bond order for a
                reaction involving these instances.

        Returns:

            A key for :attr:`_bond_orders`, which maps to the correct
            bond order.

        """

        return frozenset(map(type, functional_groups))

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        bond_orders = (
            '' if self._default_bond_orders else f'{self._bond_orders}'
        )
        return f'{self.__class__.__name__}({bond_orders})'
