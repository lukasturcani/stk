def test_caching(case_data):
    """
    Test that a :class:`.FitnessCalculator` returns cached values.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the fitness calculator to test and the
        cached fitness value.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_caching(
        fitness_calculator=case_data.fitness_calculator,
        molecule=case_data.molecule,
        fitness_value=case_data.fitness_value,
    )


def _test_caching(
    fitness_calculator,
    molecule,
    fitness_value,
):
    """
    Test that a :class:`.FitnessCalculator` returns cached values.

    Parameters
    ----------
    fitness_calculator : :class:`.FitnessCalculator`
        The fitness calculator to test.

    molecule : :class:`.Molecule`
        The molecule whose fitness value is calculated.

    fitness_value : :class:`object`
        The cached fitness value instance.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert (
        fitness_calculator.get_fitness_value(molecule) is fitness_value
    )
