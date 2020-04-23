"""
Molecule Mutation Record
========================

"""


class MutationRecord:
    """
    Abstract base class for a record of a mutation operation.

    Notes
    -----
    You might notice that the public methods of this abstract base
    class are implemented. This is just a default implementation, which
    can be used directly by users and subclasses, but can also be
    freely replaced during subclass implementation, if need be.

    """

    def __init__(self, molecule_record, mutator_name):
        """
        Initialize a :class:`.MutationRecord` instance.

        Parameters
        ----------
        molecule_record : :class:`.MoleculeRecord`
            The molecule produced by the mutation operation.

        mutator_name : :class:`str`
            The name of the mutator which carried out the mutator.

        """

        self._molecule_record = molecule_record
        self._mutator_name = mutator_name

    def get_molecule_record(self):
        """
        Get the :class:`.MoleculeRecord` produced by the mutation.

        Returns
        -------
        :class:`.MoleculeRecord`
            The molecule record.

        """

        return self._molecule_record

    def get_mutator_name(self):
        """
        Get the name of the mutator which carried out the mutation.

        Returns
        -------
        :class:`str`
            The name of the mutator.

        """

        return self._mutator_name
