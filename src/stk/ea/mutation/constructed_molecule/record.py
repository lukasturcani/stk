"""
Constructed Molecule Mutation Record
====================================

#. :class:`.ConstructedMoleculeMutationRecord`

"""

from ..molecule import MoleculeMutationRecord


class ConstructedMoleculeMutationRecord(MoleculeMutationRecord):
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
        Initialize a :class:`.ConstructedMoleculeMutationRecord`.

        Parameters
        ----------
        molecule_record : :class:`.ConstructedMoleculeRecord`
            The molecule produced by the mutation operation.

        mutator_name : :class:`str`
            The name of the mutator which carried out the mutator.

        """

        super().__init__(molecule_record, mutator_name)

    def get_molecule_record(self):
        """
        Get the molecule record produced by the mutation.

        Returns
        -------
        :class:`.ConstructedMoleculeRecord`
            The molecule record.

        """

        return super().get_molecule_record()
