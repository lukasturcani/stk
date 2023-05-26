import itertools as it

import stk


def test_select(case_data):
    """
    Test :meth:`.Selector.select`.

    This test only checks that the correct batches have been yielded,
    it does not check their order.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the selector to test and the batches which
        should be selected.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_select(
        selector=case_data.selector,
        population=case_data.population,
        selected=case_data.selected,
    )


def _test_select(selector, population, selected):
    """
    Test :meth:`.Selector.select`.

    Parameters
    ----------
    selector : :class:`.Selector`
        The selector to test.

    population : :class:`tuple` of :class:`.MoleculeRecord`
        The population from which batches are selected.

    selected : :class:`tuple` of :class:`.Batch`
            The batches which should be selected.

    Returns
    -------
    None : :class:`NoneType`

    """

    inchi = stk.Inchi()
    for batch1, batch2 in it.zip_longest(
        sorted(selector.select(population), key=get_inchis),
        sorted(selected, key=get_inchis),
    ):
        inchis1 = tuple(
            inchi.get_key(record.get_molecule()) for record in batch1
        )
        inchis2 = tuple(
            inchi.get_key(record.get_molecule()) for record in batch2
        )
        assert inchis1 == inchis2


def get_inchis(batch):
    inchi = stk.Inchi()
    return tuple(inchi.get_key(record.get_molecule()) for record in batch)
