import rdkit.Chem.AllChem as rdkit


def test_get_num_atoms(case_data):
    """
    Test :meth:`.Molecule.get_num_atoms`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the molecule to test and its correct SMILES.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_num_atoms(case_data.molecule, case_data.smiles)


def _test_get_num_atoms(molecule, smiles):
    """
    Test :meth:`.Molecule.get_num_atoms`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    smiles : :class:`str`
         The correct SMILES for `molecule`.

    Returns
    -------
    None : :class:`NoneType`

    """

    expected = rdkit.MolFromSmiles(smiles, sanitize=False)
    assert molecule.get_num_atoms() == expected.GetNumAtoms()
