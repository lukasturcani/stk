class Reaction:
    """
    An abstract base class for reactions.

    """

    def get_new_atoms(self):
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

    def get_new_bonds(self):
        """
        Yield the bonds added by the reaction.

        Yields
        ------
        :class:`.Bond`
            A bond added by the reaction.

        """

        raise NotImplementedError()

    def get_deleted_atoms(self):
        """
        Yield the atoms removed by the reaction.

        Yields
        ------
        :class:`.Atom`
            An atom deleted by the reaction.

        """

        raise NotImplementedError()
