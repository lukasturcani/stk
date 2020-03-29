"""
Molecule Crosser
================

"""


class MoleculeCrosser:
    """
    Abstract base class for molecule crossers.

    Crossers take multiple molecules and recombine them to make
    new, offspring, molecules.

    Note that despite appearances, :class:`.MoleculeCrosser` and
    :class:`.ConstructedMoleculeCrosser` are not interchangeable, you
    cannot use one where the other is required, unless explicitly
    allowed.

    Examples
    --------
    *Subclass Implementation*

    You only need to implement :meth:`._cross`. The source code of any
    of the classes listed in :mod:`.crosser` can serve as good
    examples.

    """

    def cross(self, records):
        """
        Cross `records`.

        Parameters
        ----------
        records : :class:`iterable` of :class:`.MoleculeRecord`
            The molecule records on which a crossover operation is
            performed.

        Yields
        -------
        :class:`.MoleculeCrossoverRecord`
            A record of a crossover operation.

        """

        # Can be used to decorate _cross in the future.
        yield from self._cross()

    def _cross(self, records):
        """
        Cross `records`.

        Parameters
        ----------
        molecules : :class:`iterable` of :class:`.MoleculeRecord`
            The molecule records on which a crossover operation is
            performed.

        Yields
        -------
        :class:`.MoleculeCrossoverRecord`
            A record of a crossover operation.

        """

        raise NotImplementedError()
