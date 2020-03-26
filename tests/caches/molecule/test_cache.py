from tests.utilities import is_equivalent


def test_cache(case_data):
    """
    Test a :class:`.MoleculeCache`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the cache to test and the molecule to
        place into the cache.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_cache(
        cache=case_data.cache,
        molecule=case_data.molecule,
        key=case_data.key,
    )


def _test_cache(cache, molecule, key):
    """
    Test a molecule or constructed molecule cache.

    Parameters
    ----------
    cache : :class:`.MoleculeCache`
        The cache to test.

    molecule : :class:`.Molecule`
        The molecule to put and get from the `cache`.

    key : :class:`object`
        The key used to retrieve `molecule` from the cache.

    """

    cache.put(molecule)
    retrieved = cache.get(key)
    is_equivalent(molecule, retrieved)
