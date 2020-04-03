def test_get_fitness_value(case_data):
    """
    Test :meth:`.FitnessCalculator.get_fitness_value`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the fitness calculator to test and the
        correct fitness value.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_fitness_value(
        fitness_calculator=case_data.fitness_calculator,
        molecule=case_data.molecule,
        fitness_value=case_data.fitness_value,
    )


def _test_get_fitness_value(
    fitness_calculator,
    molecule,
    fitness_value,
):
    assert (
        fitness_calculator.get_fitness_value(molecule) == fitness_value
    )
