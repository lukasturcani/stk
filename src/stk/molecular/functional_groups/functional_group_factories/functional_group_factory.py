class FunctionalGroupFactory:
    """
    Creates :class:`.FunctionalGroup` instances.

    """

    def get_functional_groups(self, molecule):
        """
        Yield functional groups in `molecule`.

        Parameters
        ----------
        :class:`.Molecule`

        Yields
        ------
        :class:`.FunctionalGroup`
            A functional group in `molecule`.

        """

        raise NotImplementedError()
