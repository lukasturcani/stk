"""
Dative Reaction Factory
=======================

"""

from .reaction_factory import ReactionFactory
from ..reactions import DativeReaction


class DativeReactionFactory(ReactionFactory):
    """
    Create reactions for :class:`.GenericFunctionalGroup` instances.

    This reaction factory assumes that the functional groups, which
    belong to the :class:`.EdgeGroup` passed to it, are
    :class:`.GenericFunctionalGroup` instances. It returns a
    :class:`.Reaction` suitable for two such instances.

    This class uses a :class:`.GenericReactionFactory` instance to
    define the :class:`Reaction`s that are performed, which are treated
    as :class:`DativeReaction`s to handle the directionality of dative
    bonds.

    Dative bonds are defined with a `bond_order` of 9.
    """

    def __init__(self, generic_reaction_factory):
        """
        Initialize a :class:`.DativeReactionFactory`.

        Parameters
        ----------
        generic_reaction_factory : :class:`.GenericReactionFactory`
            :class:`.GenericReactionFactory` that defines the
            reactions.

        """

        self._reaction_factory = generic_reaction_factory

    def get_reaction(self, construction_state, edge_group):
        reaction = self._reaction_factory.get_reaction(
            construction_state=construction_state,
            edge_group=edge_group,
        )
        return DativeReaction(reaction)

    def __str__(self):
        return repr(self)

    def __repr__(self):
        bond_orders = (
            '' if self._default_bond_orders
            else f'{self._reaction_factory._bond_orders}'
        )
        return f'{self.__class__.__name__}({bond_orders})'
