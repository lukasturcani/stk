"""
Crossover Record
=================

"""


class CrossoverRecord:
    """
    Abstract base class for a record of a crossover operation.

    Notes
    -----
    You might notice that the public methods of this abstract base
    class are implemented. This is just a default implementation, which
    can be used directly by users and subclasses, but can also be
    freely replaced during subclass implementation, if need be.

    """

    def __init__(self, molecule_record, crosser_name):
        """
        Initialize a :class:`.CrossoverRecord` instance.

        Parameters
        ----------
        molecule_record : :class:`.MoleculeRecord`
            The molecule produced by the crossover operation.

        crosser_name : :class:`str`
            The name of the crosser which carried out the crossover.

        """

        self._molecule_record = molecule_record
        self._crosser_name = crosser_name

    def get_molecule_record(self):
        """
        Get the :class:`.MoleculeRecord` produced by the crossover.

        Returns
        -------
        :class:`.MoleculeRecord`
            The molecule record.

        """

        return self._molecule_record

    def get_crosser_name(self):
        """
        Get the name of the crosser which carried out the crossover.

        Returns
        -------
        :class:`str`
            The name of the crosser.

        """

        return self._crosser_name
