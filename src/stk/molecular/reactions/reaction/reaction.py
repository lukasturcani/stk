"""
Reaction
========

.. toctree::
    :maxdepth: 2

    One One Reaction <\
stk.molecular.reactions.one_one_reaction\
>
    One Two Reaction <\
stk.molecular.reactions.one_two_reaction\
>
    Ring Amine Reaction <\
stk.molecular.reactions.ring_amine_reaction\
>
    Two Two Reaction <\
stk.molecular.reactions.two_two_reaction\
>
    Dative Reaction <\
stk.molecular.reactions.dative_reaction.dative_reaction\
>

"""

from __future__ import annotations

from collections import abc

from .new_atom import NewAtom
from .reaction_result import ReactionResult
from ...atom import Atom
from ...bond import Bond


__all__ = (
    'Reaction',
)


class Reaction:
    """
    An abstract base class for reactions.

    Reactions are used to add and remove atoms and bonds during
    :class:`.ConstructedMolecule` construction. Each subclass
    will implement a specific algorithm for doing this. Each
    :class:`.Reaction` instance should operate on a small set of
    directly relevant :class:`.FunctionalGroup` instances, and only
    modify the atoms and bonds of those instances. Normally, like 99%
    of the time, this should just be two functional groups, and you
    should ensure your topology graph only needs to react two
    functional groups at a time, if you can. However, :mod:`stk` does
    not actually care, and you can modify as many atoms and bonds as
    you want in any reaction.

    See Also:

        :mod:`.reaction_factory`
            Used for automated creation of :class:`.Reaction`
            instances. Typically, :class:`.Reaction` instances are not
            created directly, but only through some kind of
            :class:`.ReactionFactory` instance.

    Notes:

        You might notice that the public method of this abstract
        base class, :meth:`.get_result`, is implemented. This is purely
        for convenience when implementing subclasses. The implemented
        public method is simply a default implementation, which can be
        safely ignored or overridden, when implementing subclasses. The
        private methods are an implementation detail of this default
        implementation. However, they are not implemented. To use the
        default implementation of :meth:`.get_result`, the private
        methods it relies on need to be implemented.

    Examples:

        *Subclass Implementation*

        The source code of the subclasses, listed in
        :mod:`~.reaction.reaction`, can serve as good examples.

    """

    def get_result(self) -> ReactionResult:
        """
        Get the result of carrying out the reaction.

        Returns:

            Holds the results of the reaction.

        """

        return ReactionResult(
            new_atoms=tuple(self._get_new_atoms()),
            new_bonds=tuple(self._get_new_bonds()),
            deleted_atoms=tuple(self._get_deleted_atoms()),
            deleted_bonds=tuple(self._get_deleted_bonds()),
        )

    def _get_new_atoms(self) -> abc.Iterator[NewAtom]:
        """
        Yield the atoms added by the reaction.

        Yields:

            An atom added by the reaction.

        """

        raise NotImplementedError()

    def _get_new_bonds(self) -> abc.Iterator[Bond]:
        """
        Yield the bonds added by the reaction.

        Yields:

            A bond added by the reaction.

        """

        raise NotImplementedError()

    def _get_deleted_atoms(self) -> abc.Iterator[Atom]:
        """
        Yield the atoms removed by the reaction.

        Yields:

            An atom deleted by the reaction.

        """

        raise NotImplementedError()

    def _get_deleted_bonds(self) -> abc.Iterator[Bond]:
        """
        Yield the bonds removed by the reaction between existing atoms.

        Yields:

            An bond deleted by the reaction.

        """

        raise NotImplementedError()
