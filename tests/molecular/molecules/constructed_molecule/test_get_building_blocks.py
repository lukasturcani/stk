import itertools as it


def test_get_building_blocks(case_data):
    """
    Test :meth:`.ConstructedMolecule.get_building_blocks`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the constructed molecule to test and the
        correct building blocks.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_building_blocks(
        constructed_molecule=case_data.constructed_molecule,
        building_blocks=case_data.building_blocks,
    )


def _test_get_building_blocks(constructed_molecule, building_blocks):
    """
    Test :meth:`.ConstructedMolecule.get_building_blocks`.

    Parameters
    ----------
    constructed_molecule : :class:`.ConstrcutedMolecule`
        The constructed molecule to test.

    building_block : :class:`tuple` of :class:`.BuildingBlock`
        The correct building blocks.

    Returns
    -------
    None : :class:`NoneType`

    """

    for building_block1, building_block2 in it.zip_longest(
        constructed_molecule.get_building_blocks(),
        building_blocks,
    ):
        assert building_block1 is building_block2
