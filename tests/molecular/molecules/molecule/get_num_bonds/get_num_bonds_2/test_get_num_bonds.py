import rdkit.Chem.AllChem as rdkit


def test_get_num_bonds(case_data):
    """
    Test :meth:`.Molecule.get_num_bonds`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the molecule to test and its correct
        SMILES.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_num_bonds(case_data.molecule, case_data.smiles)


def _test_get_num_bonds(molecule, smiles):
    """
    Test :meth:`.Molecule.get_num_bonds`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    smiles : :class:`str`
        The correct SMILES of `molecule`.

    """

    expected = rdkit.MolFromSmiles(smiles, sanitize=False)
    assert molecule.get_num_bonds() == expected.GetNumBonds()
