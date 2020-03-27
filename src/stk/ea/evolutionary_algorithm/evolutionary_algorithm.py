"""
Evolutionary Algorithm
======================

"""


from ..fitness_normalizers import NullFitnessNormalizer
from .utilities import get_logger, get_inchi
from .implementations import Serial, Parallel


class EvolutionaryAlgorithm:
    """
    An abstract base class for evolutionary algorithms.

    Notes
    -----
    You might notice that the public methods of this abstract base
    class are implemented. This is purely for convenience, so that
    there is a default evolutionary algorithm implementation that
    users can use. However, feel free to override the default
    implementation when implementing subclasses.

    Examples
    --------
    *Subclass Implementation*

    The source code of this class can work as a good example. There
    is only one method that a subclass of
    :class:`.EvolutionaryAlgorithm` needs to implement,
    :meth:`.get_generations`, which yields :class:`.Generation`
    instances. These correspond to the generations of your
    evolutionary algorithm implementation.


    *Usage*

    There are multiple tutorials on how to use the default
    implementation of the evolutionary algorithm, which can be
    seen :ref:`here <>`

    """

    def __init__(
        self,
        initial_population,
        fitness_calculator,
        mutator,
        crosser,
        generation_selector,
        mutation_selector,
        crossover_selector,
        terminator,
        fitness_normalizer=NullFitnessNormalizer(),
        duplicate_key=get_inchi,
        logger=get_logger(),
        num_processes=1,
    ):
        """
        Initialize a :class:`EvolutionaryAlgorithm` instance.

        Parameters
        ----------

        """

        if num_processes == 1:
            self._implementation = Serial(
                initial_population=initial_population,
                fitness_calculator=fitness_calculator,
                mutator=mutator,
                crosser=crosser,
                generation_selector=generation_selector,
                mutation_selector=mutation_selector,
                crossover_selector=crossover_selector,
                terminator=terminator,
                fitness_normalizer=fitness_normalizer,
                duplicate_key=duplicate_key,
                logger=logger,
            )

        else:
            self._implementation = Parallel(
                initial_population=initial_population,
                fitness_calculator=fitness_calculator,
                mutator=mutator,
                crosser=crosser,
                generation_selector=generation_selector,
                mutation_selector=mutation_selector,
                crossover_selector=crossover_selector,
                terminator=terminator,
                fitness_normalizer=fitness_normalizer,
                duplicate_key=duplicate_key,
                logger=logger,
                num_processes=num_processes,
            )

    def get_generations(self):
        """
        Yield the generations of the evolutionary algorithm.

        Yields
        ------
        :class:`.Generation`
            A generation.

        """

        yield from self._implementation.get_generations(self)
