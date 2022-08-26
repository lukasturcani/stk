def test_get_bond_infos(case_data):
    """
    Test :meth:`.ConstructedMolecule.get_bond_infos`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the constructed molecule to test and the
        correct number of new bonds.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_bond_infos(
        constructed_molecule=case_data.constructed_molecule,
        num_new_bonds=case_data.num_new_bonds,
    )


def _test_get_bond_infos(constructed_molecule, num_new_bonds):
    """
    Test :meth:`.ConstructedMolecule.get_bond_infos`.

    Parameters
    ----------
    constructed_molecule : :class:`.ConstructedMolecule`
        The constructed molecule to test.

    num_new_bonds : :class:`int`
        The correct number of new bonds added.

    Returns
    -------
    None : :class:`NoneType`


    """

    new_bonds = filter(
        lambda bond_info: bond_info.get_building_block() is None,
        constructed_molecule.get_bond_infos(),
    )
    assert sum(1 for _ in new_bonds) == num_new_bonds
    assert constructed_molecule.get_num_bonds() == sum(
        1 for _ in constructed_molecule.get_bond_infos()
    )
