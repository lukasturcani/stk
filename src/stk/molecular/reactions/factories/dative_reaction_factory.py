"""
Dative Reaction Factory
=======================

"""

from ..reactions import DativeReaction
from .reaction_factory import ReactionFactory


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

    def __init__(self, reaction_factory):
        """
        Initialize a :class:`.DativeReactionFactory`.

        Parameters
        ----------
        reaction_factory : :class:`.GenericReactionFactory`
            Used to create reactions.

        """

        self._reaction_factory = reaction_factory

    def get_reaction(self, construction_state, edge_group):
        reaction = self._reaction_factory.get_reaction(
            construction_state=construction_state,
            edge_group=edge_group,
        )
        return DativeReaction(reaction)

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return (
            f'{self.__class__.__name__}'
            f'({self._reaction_factory})'
        )
