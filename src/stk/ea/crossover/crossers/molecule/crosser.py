"""
Molecule Crosser
================

.. toctree::
    :maxdepth: 2

    Genetic Recombination <\
stk.ea.crossover.crossers.molecule.genetic_recombination\
>

"""


class MoleculeCrosser:
    """
    Abstract base class for molecule crossers.

    Crossers take multiple molecules and recombine them to make
    new, offspring, molecules.

    Examples
    --------
    *Subclass Implementation*

    You only need to implement :meth:`.cross`. The source code of any
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
        :class:`.CrossoverRecord`
            A record of a crossover operation.

        """

        raise NotImplementedError()
