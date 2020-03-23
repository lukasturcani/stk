def test_to_json(case_data):
    """
    Test :meth:`.to_json`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the jsonizer and the correct json of a
        molecule.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_to_json(
        jsonizer=case_data.jsonizer,
        molecule=case_data.molecule,
        json=case_data.json,
    )


def _test_to_json(jsonizer, molecule, json):
    """
    Test :meth:`.to_json`.

    Parameters
    ----------
    jsonizer : :class:`.MoleculeJsonizer` or \
            :class:`.ConstructedMoleculeJsonizer`
        The jsonizer to test.

    molecule : :class:`.Molecule`
        The molecule to JSONize.

    json : :class:`dict`
        The correct JSON of `molecule`.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert jsonizer.to_json(molecule) == json
