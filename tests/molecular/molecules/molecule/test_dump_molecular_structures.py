def test_dump_molecular_structures(
    request,
    tmp_path,
    molecule_db,
    name_db,
    case_data,
):
    """
    Dump molecular structures.

    This test dumps molecules to files and to a MongoDB database
    so that they can be visually inspected.

    Parameters
    ----------
    request : :class:`pytest.FixtureRequest`
        Holds information about the requesting test.

    tmp_path : :class:`pathlib2.Path`
        A path into which the structure of the molecule is saved.

    molecule_db : :class:`.MoleculeDatabase`
        A database into which the structure of the molecule  is saved.

    name_db : :class:`.ValueDatabase`
        A database into which name the name of the molecule is saved.

    Returns
    -------
    None : :class:`NoneType`

    """

    case_data.molecule.write(tmp_path / f"{request.node.name}.mol")
    molecule_db.put(case_data.molecule)
    name_db.put(case_data.molecule, request.node.name)
