from typing import NamedTuple


class ReactionResult(NamedTuple):
    """
    The result of a reaction.

    Attributes
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

    new_atoms: object
    new_bonds: object
    deleted_atoms: object


class Reaction:
    """
    An abstract base class for reactions.

    """

    def get_result(self):
        """
        Get the result of carrying out the reaction.

        Returns
        -------
        :class:`._ReactionResult`
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
