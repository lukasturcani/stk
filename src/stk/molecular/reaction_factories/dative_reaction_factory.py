"""
Dative Reaction Factory
=======================

"""

from .reaction_factory import ReactionFactory
from .generic_reaction_factory import GenericReactionFactory
from ..construction_state import ConstructionState
from ..edge_group import EdgeGroup
from ..reactions import DativeReaction


__all__ = (
    'DativeReactionFactory',
)


class DativeReactionFactory(ReactionFactory):
    """
    Create :class:`.DativeReaction` instances.

    This reaction factory assumes that the functional groups, which
    belong to the :class:`.EdgeGroup` passed to it, are
    :class:`.GenericFunctionalGroup` instances. It returns a
    :class:`.Reaction` suitable for two such instances.

    Dative bonds are defined with a `bond_order` of 9,
    running from the non-metal atom to the metal atom.

    """

    def __init__(
        self,
        reaction_factory: GenericReactionFactory,
    ) -> None:
        """
        Initialize a :class:`.DativeReactionFactory`.

        Parameters:

            reaction_factory:
                Used to create reactions.

        """

        self._reaction_factory = reaction_factory

    def get_reaction(
        self,
        construction_state: ConstructionState,
        edge_group: EdgeGroup,
    ) -> DativeReaction:

        reaction = self._reaction_factory.get_reaction(
            construction_state=construction_state,
            edge_group=edge_group,
        )
        return DativeReaction(reaction)

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return (
            f'{self.__class__.__name__}'
            f'({self._reaction_factory})'
        )
