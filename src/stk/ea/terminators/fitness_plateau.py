
class FitnessPlateau(Terminator):
    """
    Checks if the fittest molecules remain the same.

    """

    def __init__(self, num_generations, top_members=1):
        """
        Initialize a :class:`FitnessPlateau` instance.

        Parameters
        ----------
        num_generations : :class:`int`
            Number of generations in which the fittest molecules did
            not change.

        top_members : :class:`int`, optional
            The number of fittest molecules which are checked. This
            number needs to be smaller than the population size.

        """

        self._num_generations = num_generations
        self._top_members = top_members

    def terminate(self, progress):
        """
        Check if the fittest molecules changed between generations.

        Parameters
        ----------
        progress : :class:`.Population`
            A population where every generation is a subpopulation.

        Returns
        -------
        :class:`bool`
            ``True`` if the :attr:`top_members` have not changed over
            the last :attr:`num_gens` generations.

        """

        # Check that the EA has run for more than num_gens generations.
        if len(progress.subpopulations) >= self._num_generations:
            gens = set()
            for i in range(self._num_generations):
                gen = sorted(
                    progress.subpopulations[-i-1],
                    reverse=True,
                    key=lambda mol: mol.fitness
                )
                # Get the top members of the generation.
                keys = frozenset(
                    mol.get_identity_key() for mol in gen[:self._top_members]
                )
                gens.add(keys)
            unique_gens = len(gens)
            if unique_gens == 1:
                return True
            else:
                return False

        return False
