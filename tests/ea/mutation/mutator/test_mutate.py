from tests.utilities import is_equivalent


def test_mutate(case_data):
    """
    Test :meth:`mutate`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the mutator to test and the correct mutation
        record.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_mutate(
        mutator=case_data.mutator,
        record=case_data.record,
        mutation_record=case_data.mutation_record,
    )


def _test_mutate(mutator, record, mutation_record):
    """
    Test :meth:`mutate`.

    Parameters
    ----------
    mutator : :class:`.MoleculeMutator`
        The mutator to test.

    record : :class:`.MoleculeRecord`
        The molecule to mutate.

    mutation_record : :class:`.MutationRecord`
        The correct mutation record.

    Returns
    -------
    None : :class:`NoneType`

    """

    result = mutator.mutate(record)
    assert result.get_mutator_name() == mutation_record.get_mutator_name()
    is_equivalent(
        result.get_molecule_record().get_molecule(),
        mutation_record.get_molecule_record().get_molecule(),
    )
