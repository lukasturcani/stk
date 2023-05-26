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
    """
    Test :meth:`.FitnessCalculator.get_fitness_value`.

    Parameters
    ----------
    fitness_calculator : :class:`.FitnessCalculator`
        The fitness calculator to test.

    molecule : :class:`.Molecule`
        The molecule whose fitness is calculated.

    fitness_value : :class:`object`
        The correct fitness value.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert fitness_calculator.get_fitness_value(molecule) == fitness_value
