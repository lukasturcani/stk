"""
Generation
==========

"""


class Generation:
    """
    A generation of the :class:`.EvolutionaryAlgorithm`.

    """

    def __init__(
        self,
        id,
        molecule_records,
        mutation_records,
        crossover_records,
    ):
        self._id = id
        self._molecule_records = molecule_records
        self._mutation_records = mutation_records
        self._crossover_records = crossover_records

    def get_id(self):
        """
        Get the id of the generation.

        Returns
        -------
        :class:`int`
            The id of the generation.

        """

        return self._id

    def get_molecule_records(self):
        """
        Yield the molecule records in the generation.

        Yields
        ------
        :class:`.MoleculeRecord`
            A molecule record.

        """

        yield from self._molecule_records

    def get_mutation_records(self):
        """
        Yield the mutation records in the generation.

        Yields
        ------
        :class:`.MutationRecord`
            A mutation record.

        """

        yield from self._mutation_records

    def get_crossover_records(self):
        """
        Yield the crossover records in the generation.

        Yields
        ------
        :class:`.CrossoverRecord`
            A crossover record.

        """

        yield from self._crossover_records
