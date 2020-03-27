"""
Molecule Crosser
================

"""


class MoleculeCrosser:
    """
    Abstract base class for molecule crossers.

    Crossers take multiple molecules and recombine them to make
    new, offspring, molecules.

    """

    def cross(self, molecules):
        """
        Cross `molecules`.

        Parameters
        ----------
        molecules : :class:`iterable` of :class:`.MoleculeRecord`
            The molecules on which a crossover operation is performed.

        Yields
        -------
        :class:`.CrossoverRecord`
            A record of a crossover operation.

        """

        # Can be used to decorate _cross in the future.
        yield from self._cross()

    def _cross(self, molecules):
        """
        Cross `molecules`.

        Parameters
        ----------
        molecules : :class:`iterable` of :class:`.MoleculeRecord`
            The molecules on which a crossover operation is performed.

        Yields
        -------
        :class:`.CrossoverRecord`
            A record of a crossover operation.

        """

        raise NotImplementedError()
