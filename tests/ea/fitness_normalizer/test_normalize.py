import itertools as it
import numpy as np

from tests.utilities import is_equivalent


def test_normalize(case_data):
    """
    Test :meth:`.FitnessNormalizer.normalize`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the fitness normalizer to test and the
        normalized population.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_normalize(
        fitness_normalizer=case_data.fitness_normalizer,
        population=case_data.population,
        normalized=case_data.normalized,
    )


def _test_normalize(
    fitness_normalizer,
    population,
    normalized,
):
    """
    Test :meth:`.FitnessNormalizer.normalize`.

    Parameters
    ----------
    fitness_normalizer : :class:`.FitnessNormalizer`
        The fitness normalizer to test.

    population : :class:`tuple` of :class:`.MoleculeRecord`
        The population which is normalized.

    normalized : :class:`tuple` of :class:`.MoleculeRecord`
        The normalized `population`.

    Returns
    -------
    None : :class:`NoneType`

    """

    for record1, record2 in it.zip_longest(
        fitness_normalizer.normalize(population),
        normalized,
    ):
        is_equivalent(
            record1.get_molecule(),
            record2.get_molecule(),
        )
        fitness_value = record1.get_fitness_value()
        if isinstance(fitness_value, np.ndarray):
            fitness_value = tuple(fitness_value)

        assert (fitness_value == record2.get_fitness_value())
        assert (
            record1.get_fitness_value(normalized=False)
            == record2.get_fitness_value(normalized=False)
        )
