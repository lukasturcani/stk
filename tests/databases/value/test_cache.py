def test_cache(case_data):
    """
    Test a cache.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the cache to test and the value to put
        into it.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_cache(
        cache=case_data.cache,
        molecule=case_data.molecule,
        value=case_data.value,
    )


def _test_cache(cache, molecule, value):
    """
    Test a cache.

    Parameters
    ----------
    cache : class:`.MoleculeValueCache` or \
            :class:`.ConstructedMoleculeValueCache`
        The cache to test.

    molecule : :class:`.Molecule`
        The molecule to test.

    value : :class:`object`
        The value to put into the cache.

    Returns
    -------
    None : :class:`NoneType`

    """

    cache.put(molecule, value)
    assert cache.get(molecule) == value
