"""
Constructed Molecule Mutator
============================

#. :class:`.RandomBuildingBlock`
#. :class:`.RandomTopologyGraph`
#. :class:`.SimilarBuildingBlock`

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

    You only need to implement :meth:`.mutate`. The source code of any
    of the classes listed in :mod:`.mutator` can serve as good
    examples.

    """

    def mutate(self, record):
        """
        Return a mutant of `record`.

        Parameters
        ----------
        record : :class:`.ConstructedMoleculeRecord`
            The molecule to be mutated.

        Returns
        -------
        :class:`.ConstructedMoleculeMutationRecord`
            A record of the mutation.

        """

        raise NotImplementedError()
