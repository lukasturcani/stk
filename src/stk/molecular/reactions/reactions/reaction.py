class Reaction:
    """
    An abstract base class for reactions.

    """

    def get_new_atoms(self):
        """
        Yield the atoms added by the reaction.

        Yields
        ------
        :class:`.Atom`
            An atom added by the reaction.

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
