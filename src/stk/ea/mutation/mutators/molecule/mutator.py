"""
Molecule Mutator
================

.. toctree::
    :maxdepth: 2

    Random Building Block <\
stk.ea.mutation.mutators.molecule.random_building_block\
>
    Random Topology Graph <\
stk.ea.mutation.mutators.molecule.random_topology_graph\
>
    Similar Building Block <\
stk.ea.mutation.mutators.molecule.similar_building_block\
>

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

        Raises
        ------
        :class:`.MutationPreconditionViolation`
            If the molecule which is meant to be mutated cannot be,
            because it does not satisfy the necessary preconditions for
            the mutation operation. Reading the documentation of
            :class:`.MutationPreconditionViolation`
            is strongly strongly recommended for understanding why
            you might want to raise this error.

        """

        raise NotImplementedError()
