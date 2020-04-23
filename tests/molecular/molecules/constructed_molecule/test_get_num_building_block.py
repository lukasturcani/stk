def test_get_num_building_block(case_data):
    """
    Test :meth:`.ConstructedMolecule.get_num_building_block`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the constructed molecule to test and the
        correct number of each building block.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_num_building_block(
        constructed_molecule=case_data.constructed_molecule,
        num_building_blocks=case_data.num_building_blocks,
    )


def _test_get_num_building_block(
    constructed_molecule,
    num_building_blocks,
):
    """
    Test :meth:`.ConstructedMolecule.get_num_building_blocks`.

    Parameters
    ----------
    constructed_molecule : :class:`.ConstructedMolecule`
        The constructed molecule to test.

    num_building_blocks : :class:`dict`
        Maps each building block in `constructed_molecule` to the
        number of times it was used during construction.

    Returns
    -------
    None : :class:`NoneType`

    """

    for building_block in constructed_molecule.get_building_blocks():
        assert (
            constructed_molecule.get_num_building_block(building_block)
            == num_building_blocks[building_block]
        )
