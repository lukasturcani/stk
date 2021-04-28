import rdkit.Chem.AllChem as rdkit


def test_get_atoms(
    request,
    tmp_path,
    molecule_db,
    name_db,
    case_data,
):
    """
    Test :meth:`.Molecule.get_atoms`.

    Parameters
    ----------
    request : :class:`pytest.FixtureRequest`
        Holds information about the requesting test function.

    tmp_path : :class:`pathlib2.Path`
        A path into which the structure of the molecule being tested
        is saved.

    molecule_db : :class:`.MoleculeDatabase`
        A database into which the molecule being tested is saved.

    name_db : :class:`.ValueDatabase`
        A database into which name the molecule being tested is
        saved.

    case_data : :class:`.CaseData`
        A test case. Holds the molecule to test and the correct
        SMILES.

    Returns
    -------
    None : :class:`NoneType`

    """

    # Visual inspection of the molecule can be useful.
    case_data.molecule.write(tmp_path / 'molecule.mol')
    molecule_db.put(case_data.molecule)
    name_db.put(case_data.molecule, request.node.name)
    _test_get_atoms(case_data.molecule, case_data.smiles)


def _test_get_atoms(molecule, smiles):
    """
    Test :meth:`.Molecule.get_atoms`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    smiles : :class:`str`
        The correct SMILES of `molecule`.

    Returns
    -------
    None : :class:`NoneType`

    """

    # get_rdkit_mol() tests get_atoms() and get_bonds(), if those are
    # broken, the produced rdkit molecule will not have the expected
    # smiles.
    result = rdkit.MolToSmiles(get_rdkit_mol(molecule))
    # Printing is useful for debugging.
    print(result)
    assert result == smiles


def get_rdkit_mol(molecule):
    rdkit_mol = rdkit.EditableMol(rdkit.Mol())
    for atom in molecule.get_atoms():
        rdkit_atom = rdkit.Atom(atom.get_atomic_number())
        rdkit_atom.SetFormalCharge(atom.get_charge())
        rdkit_mol.AddAtom(rdkit_atom)

    for bond in molecule.get_bonds():
        rdkit_mol.AddBond(
            beginAtomIdx=bond.get_atom1().get_id(),
            endAtomIdx=bond.get_atom2().get_id(),
            order=(
                rdkit.BondType.DATIVE if bond.get_order() == 9
                else rdkit.BondType(bond.get_order())
            ),
        )
    return rdkit_mol.GetMol()
