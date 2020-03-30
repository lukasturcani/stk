"""
Molecule Mutator
================

"""


class ConstructedMoleculeMutator:
    """
    Abstract base class for molecule mutators.

    Note that despite appearances, :class:`.MoleculeMutator` and
    :class:`.ConstructedMoleculeMutator` are not interchangeable, you
    cannot use one where the other is required, unless explicitly
    allowed.

    Examples
    --------
    *Subclass Implementation*

    You only need to implement :meth:`._mutate`. The source code of any
    of the classes listed in :mod:`.mutator` can serve as good
    examples.

    """

    def mutate(self, record):
        """
        Return a mutant of `record`.

        Parameters
        ----------
        record : :class:`.MoleculeRecord`
            The molecule to be mutated.

        Returns
        -------
        :class:`.MoleculeMutationRecord`
            A record of the mutation.

        """

        # Can be used to decorate _mutate in the future.
        return self._mutate(record)

    def _mutate(self, mol):
        """
        Return a mutant of `record`.

        Parameters
        ----------
        record : :class:`.MoleculeRecord`
            The molecule to be mutated.

        Returns
        -------
        :class:`.MoleculeMutationRecord`
            A record of the mutation.

        """

        raise NotImplementedError()
