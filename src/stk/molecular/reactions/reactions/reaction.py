"""
Reaction
========

#. :class:`.OneOneReaction`
#. :class:`.OneTwoReaction`
#. :class:`.RingAmineReaction`
#. :class:`.TwoTwoReaction`

"""


class ReactionResult:
    """
    The result of a reaction.

    """

    __slots__ = ['_new_atoms', '_new_bonds', '_deleted_atoms']

    def __init__(self, new_atoms, new_bonds, deleted_atoms):
        """
        Initialize a :class:`.ReactionResult` instance.

        Parameters
        ----------
        new_atoms : :class:`tuple`
            A :class:`tuple` of the form ``(atom, position)``, holding
            an :class:`.Atom` added
            by the reaction and its position as a
            :class:`numpy.ndarray`. New atoms have a negative id, and
            will need to be assigned a new one when added to the
            :class:`.ConstructedMolecule`.

        new_bonds : :class:`tuple` of :class:`.Bond`
            The bonds added by the reaction.

        deleted_atoms : :class:`tuple` of :class:`.Atom`
            The atoms deleted by the reaction.

        """

        self._new_atoms = new_atoms
        self._new_bonds = new_bonds
        self._deleted_atoms = deleted_atoms

    def get_new_atoms(self):
        """
        Get the new atoms added by the reaction.

        Returns
        -------
        :class:`tuple`
            A :class:`tuple` of the form ``(atom, position)``, holding
            an :class:`.Atom` added
            by the reaction and its position as a
            :class:`numpy.ndarray`. New atoms have a negative id, and
            will need to be assigned a new one when added to the
            :class:`.ConstructedMolecule`.

        """

        return self._new_atoms

    def get_new_bonds(self):
        """
        Get the new bonds added by the reaction.

        Returns
        -------
        :class:`tuple` of :class:`.Bond`
            The new bonds added by the reaction.

        """
        return self._new_bonds

    def get_deleted_atoms(self):
        """
        Get the atoms deleted by the reaction.

        Returns
        -------
        :class:`tuple` of :class:`.Atom`
            The atoms deleted by the reaction.

        """

        return self._deleted_atoms


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

    See Also
    --------
    :mod:`.reaction_factory`
        Used for automated creation of :class:`.Reaction` instances.
        Typically, :class:`.Reaction` instances are not created
        directly, but only through some kind of
        :class:`.ReactionFactory` instance.

    Notes
    -----
    You might notice that the public method of this abstract
    base class, :meth:`.get_result`, is implemented. This is purely for
    convenience when implementing subclasses. The implemented public
    method is simply a default implementation, which can be safely
    ignored or overridden, when implementing subclasses. The private
    methods are an implementation detail of this default
    implementation. However, they are not implemented. To use the
    default implementation of :meth:`.get_result`, the private methods
    it relies on need to be implemented.

    Examples
    --------
    *Subclass Implementation*

    The source code of the subclasses, listed in :mod:`.reaction`, can
    serve as good examples.

    """

    def get_result(self):
        """
        Get the result of carrying out the reaction.

        Returns
        -------
        :class:`.ReactionResult`
            Holds the results of the reaction.

        """

        return ReactionResult(
            new_atoms=tuple(self._get_new_atoms()),
            new_bonds=tuple(self._get_new_bonds()),
            deleted_atoms=tuple(self._get_deleted_atoms()),
        )

    def _get_new_atoms(self):
        """
        Yield the atoms added by the reaction.

        Yields
        ------
        :class:`tuple`
            A :class:`tuple` of the form ``(atom, position)``, holding
            an :class:`.Atom` added
            by the reaction and its position as a
            :class:`numpy.ndarray`. New atoms have a negative id, and
            will need to be assigned a new one when added to the
            :class:`.ConstructedMolecule`.

        """

        raise NotImplementedError()

    def _get_new_bonds(self):
        """
        Yield the bonds added by the reaction.

        Yields
        ------
        :class:`.Bond`
            A bond added by the reaction.

        """

        raise NotImplementedError()

    def _get_deleted_atoms(self):
        """
        Yield the atoms removed by the reaction.

        Yields
        ------
        :class:`.Atom`
            An atom deleted by the reaction.

        """

        raise NotImplementedError()
