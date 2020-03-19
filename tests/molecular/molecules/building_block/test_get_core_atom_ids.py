import itertools as it


def test_get_core_atom_ids(case_data):
    """
    Test :meth:`.BuildingBlock.get_core_atom_ids`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the building block to test and the correct
        core atom ids.

    Returns
    -------
    None : :class:`NoneType`

    """

    for id1, id2 in it.zip_longest(
        case_data.building_block.get_core_atom_ids(),
        case_data.core_atom_ids,
    ):
        assert id1 == id2
