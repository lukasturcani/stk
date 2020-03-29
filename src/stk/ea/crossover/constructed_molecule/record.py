"""
Constructed Molecule Crossover Record
=====================================

#. :class:`.ConstructedMoleculeCrossoverRecord`

"""


from ..molecule import MoleculeCrossoverRecord


class ConstructedMoleculeCrossoverRecord(MoleculeCrossoverRecord):
    """
    Abstract base class for a record of a crossover operation.

    Notes
    -----
    You might notice that the public methods of this abstract base
    class are implemented. This is just a default implementation, which
    can be used directly by users and subclasses, but can also be
    freely replaced during subclass implementation, if need be.

    """

    def __init__(
        self,
        molecule_record,
        crosser_name,
    ):
        """
        Initialize a :class:`.ConstructedMoleculeCrossoverRecord`.

        Parameters
        ----------
        molecule_record : :class:`.ConstructedMoleculeRecord`
            The molecule produced by the crossover operation.

        crosser_name : :class:`str`
            The name of the crosser which carried out the crossover.

        """

        super().__init__(molecule_record, crosser_name)

    def get_molecule_record(self):
        """
        Get the molecule record produced by the crossover.

        Returns
        -------
        :class:`.ConstructedMoleculeRecord`
            The molecule record.

        """

        yield from self._molecule_record
