def test_to_rdkit_mol(molecule):
    """
    Test :meth:`.Molecule.to_rdkit_mol`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    Returns
    -------
    None : :class:`NoneType`

    """

    rdkit_molecule = molecule.to_rdkit_mol()
    assert rdkit_molecule.GetNumConformers() == 1

    assert rdkit_molecule.GetNumAtoms() == molecule.get_num_atoms()
    atoms = zip(molecule.get_atoms(), rdkit_molecule.GetAtoms())
    for atom, rdkit_atom in atoms:
        assert atom.get_charge() == rdkit_atom.GetFormalCharge()
        assert atom.get_atomic_number() == rdkit_atom.GetAtomicNum()

    assert molecule.get_num_bonds() == rdkit_molecule.GetNumBonds()
    bonds = zip(molecule.get_bonds(), rdkit_molecule.GetBonds())
    for bond, rdkit_bond in bonds:
        assert bond.get_order() == rdkit_bond.GetBondTypeAsDouble()
        assert bond.get_atom1().get_id() == rdkit_bond.GetBeginAtomIdx()
        assert bond.get_atom2().get_id() == rdkit_bond.GetEndAtomIdx()
