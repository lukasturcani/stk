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

        None : :class:`NoneType`
            If `record` cannot be mutated.

        """

        raise NotImplementedError()
