def test_get_atom_infos(case_data):
    """
    Test :meth:`.ConstructedMolecule.get_atom_infos`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the constructed molecule to test and the
        correct number of new atoms.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_atom_infos(
        constructed_molecule=case_data.constructed_molecule,
        num_new_atoms=case_data.num_new_atoms,
    )


def _test_get_atom_infos(constructed_molecule, num_new_atoms):
    """
    Test :meth:`.ConstructedMolecule.get_atom_infos`.

    Parameters
    ----------
    constructed_molecule : :class:`.ConstructedMolecule`
        The constructed molecule to test.

    num_new_atoms : :class:`int`
        The correct number of new atoms.

    Returns
    -------
    None : :class:`NoneType`

    """

    new_atoms = filter(
        lambda atom_info: atom_info.get_building_block() is None,
        constructed_molecule.get_atom_infos(),
    )
    assert sum(1 for _ in new_atoms) == num_new_atoms
    assert constructed_molecule.get_num_atoms() == sum(
        1 for _ in constructed_molecule.get_atom_infos()
    )
