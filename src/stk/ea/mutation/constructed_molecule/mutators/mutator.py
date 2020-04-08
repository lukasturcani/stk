"""
Constructed Molecule Mutator
============================

#. :class:`.RandomBuildingBlock`
#. :class:`.RandomTopologyGraph`
#. :class:`.SimilarBuildingBlock`

"""


class MoleculeMutator:
    """
    Abstract base class for molecule mutators.

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
        record : :class:`.MoleculeRecord`
            The molecule to be mutated.

        Returns
        -------
        :class:`.MutationRecord`
            A record of the mutation.

        """

        raise NotImplementedError()
