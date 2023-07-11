from typing import Any

import numpy as np
import stk

from .case_data import CaseData


def test_normalize(case_data: CaseData) -> None:
    _test_normalize(
        fitness_normalizer=case_data.fitness_normalizer,
        population=case_data.fitness_values,
        normalized=case_data.normalized,
    )


def _test_normalize(
    fitness_normalizer: stk.FitnessNormalizer[stk.MoleculeRecord[Any]],
    population: dict[stk.MoleculeRecord[Any], Any],
    normalized: dict[stk.MoleculeRecord[Any], Any],
) -> None:
    result = fitness_normalizer.normalize(population)
    assert len(result) == len(population)
    for record, fitness_value in result.items():
        if isinstance(fitness_value, np.ndarray):
            fitness_value = tuple(fitness_value)
        assert fitness_value == normalized[record]
