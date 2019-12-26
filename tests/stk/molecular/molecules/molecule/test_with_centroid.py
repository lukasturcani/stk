def test_with_centroid(molecule, get_atom_ids):
    # Keep clone to test for immutability.
    clone = molecule.clone()

    new = molecule.with_centroid(
        position=[1, 2, 3],
        atom_ids=get_atom_ids(molecule),
    )
    assert np.allclose(
        a=new.get_centroid(atom_ids=get_atom_ids(molecule)),
        b=[1, 2, 3],
        atol=1e-32,
    )
    _test_unchanged(clone, molecule)
