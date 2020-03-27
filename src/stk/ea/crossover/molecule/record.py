"""
Crossover Record
================

"""


class CrossoverRecord:
    """
    Holds a record of a crossover operation.

    """

    def __init__(self, molecule_record, crosser_name):
        self._molecule_record = molecule_record
        self._crosser_name = crosser_name

    def get_molecule_record(self):
        """
        Get the :class:`.MoleculeRecord` produced by the operation.

        Returns
        -------
        :class:`.MoleculeRecord`
            The molecule record.

        """

        yield from self._molecule_record

    def get_crosser_name(self):
        """
        Get the name of the crosser which created this record.

        Returns
        -------
        :class:`str`
            The name of the crosser.

        """

        return self._crosser_name
