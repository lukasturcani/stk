import typing

import stk

from .case_data import CaseData


def test_caching(case_data: CaseData) -> None:
    """
    Test that a :class:`.FitnessCalculator` returns cached values.
    """

    _test_caching(
        fitness_calculator=case_data.fitness_calculator,
        record=case_data.record,
        fitness_value=case_data.fitness_value,
    )


def _test_caching(
    fitness_calculator: stk.FitnessCalculator,
    record: stk.MoleculeRecord,
    fitness_value: typing.Any,
):
    assert fitness_calculator.get_fitness_value(record) is fitness_value
