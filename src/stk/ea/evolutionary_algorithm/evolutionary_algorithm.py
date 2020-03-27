"""
Evolutionary Algorithm
======================

"""


import logging
import os

from ..fitness_normalizers import NullFitnessNormalizer
from .utilities import get_logger




class EvolutionaryAlgorithm:
    """
    An abstract base class for evolutionary algorithms.

    Notes
    -----


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
        logger=get_logger(),
    ):
        """
        Initialize a :class:`EvolutionaryAlgorithm` instance.

        Parameters
        ----------

        """

        self._population = initial_population
        self._fitness_calculator = fitness_calculator
        self._fitness_normalizer = fitness_normalizer
        self._mutator = mutator
        self._crosser = crosser
        self._generation_selector = generation_selector
        self._mutation_selector = mutation_selector
        self._crossover_selector = crossover_selector
        self._terminator = terminator

    def get_generations(self):
        """
        Yield the generations of the evolutionary algorithm.

        Yields
        ------
        :class:`.Generation`
            A generation.

        """




# Define the formatter for logging messages.
try:
    f = '\n' + '='*os.get_terminal_size().columns + '\n\n'
except OSError:
    # When testing os.get_terminal_size() will fail because stdout is
    # not connected to a terminal.
    f = '\n' + '='*100 + '\n\n'
formatter = logging.Formatter(
    fmt=f+('%(asctime)s - %(levelname)s - %(name)s - %(message)s'),
    datefmt='%H:%M:%S'
)

# Define logging handlers.
errorhandler = logging.FileHandler(
    'output/scratch/errors.log',
    delay=True
)
errorhandler.setLevel(logging.ERROR)

streamhandler = logging.StreamHandler()

errorhandler.setFormatter(formatter)
streamhandler.setFormatter(formatter)

# Get the loggers.
rootlogger = logging.getLogger()
rootlogger.addHandler(errorhandler)
rootlogger.addHandler(streamhandler)


logger = logging.getLogger(__name__)
