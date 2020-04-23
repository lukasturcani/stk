def test_get_key(case_data):
    """
    Test :meth:`.get_key`.

    Parameters
    ----------
    case_data : :class:`.CaseDat`
        A test case. Holds the key maker to test and the correct key
        it should produce.

    Returns
    -------
    None : :class`NoneType`

    """

    _test_get_key(
        key_maker=case_data.key_maker,
        molecule=case_data.molecule,
        key=case_data.key,
    )


def _test_get_key(key_maker, molecule, key):
    """
    Test :meth:`.get_key`.

    Parameters
    ----------
    key_maker : :class:`.MoleculeKeyMaker`
        The key maker to test.

    molecule : :class:`.Molecule`
        The molecule to pass to the `key_maker`.

    key : :class:`object`
        The correct key of `molecule`.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert key_maker.get_key(molecule) == key
