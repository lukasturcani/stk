
class Mutator(EAOperation):
    """
    Creates mutants.

    """

    def mutate(self, mol):
        """
        Return a mutant of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be mutated.

        Returns
        -------
        mol : :class:`.Molecule`
            The mutant.

        """

        # Can be used to decorate _mutate in the future.
        return self._mutate(mol)

    def _mutate(self, mol):
        """
        Return a mutant of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be mutated.

        Returns
        -------
        mol : :class:`.Molecule`
            The mutant.

        Raises
        ------
        :class:`NotImplementedError`
            This is a virtual method which must be implemented by
            a subclass.

        """

        raise NotImplementedError()
