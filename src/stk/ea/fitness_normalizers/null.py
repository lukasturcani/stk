
class NullFitnessNormalizer(FitnessNormalizer):
    """
    Does nothing.

    """

    def _normalize(self, population):
        return population.get_fitness_values()


