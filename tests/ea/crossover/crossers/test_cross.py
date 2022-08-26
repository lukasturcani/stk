from tests.utilities import is_equivalent


def test_cross(case_data):
    """
    Test :meth:`cross`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the crosser to test and the correct
        crossover records.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_cross(
        crosser=case_data.crosser,
        records=case_data.records,
        crossover_records=case_data.crossover_records,
    )


def _test_cross(crosser, records, crossover_records):
    """
    Test :meth:`cross`.

    Parameters
    ----------
    crosser : :class:`.MoleculeCrosser`
        The crosser to test.

    records : :class:`tuple` of :class:`.MoleculeRecord`
        The molecules to cross.

    crossover_records : :class:`tuple` of :class:`.CrossoverRecord`
        The correct offspring.

    Returns
    -------
    None : :class:`NoneType`

    """

    for record1, record2 in zip(
        crosser.cross(records),
        crossover_records,
    ):
        assert record1.get_crosser_name() == record2.get_crosser_name()
        is_equivalent(
            record1.get_molecule_record().get_molecule(),
            record2.get_molecule_record().get_molecule(),
        )
