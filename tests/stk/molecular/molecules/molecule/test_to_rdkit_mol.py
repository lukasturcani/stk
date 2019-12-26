def test_to_rdkit_mol(molecule):
    rdkit_molecule = molecule.to_rdkit_mol()
    assert rdkit_molecule.GetNumConformers() == 1

    assert rdkit_molecule.GetNumAtoms() == molecule.get_num_atoms()
    atoms = zip(molecule.get_atoms(), rdkit_molecule.GetAtoms())
    for atom, rdkit_atom in atoms:
        assert atom.charge == rdkit_atom.GetFormalCharge()
        assert atom.atomic_number == rdkit_atom.GetAtomicNum()
        assert atom.mass == rdkit_atom.GetMass()

    assert molecule.get_num_bonds() == rdkit_molecule.GetNumBonds()
    bonds = zip(molecule.get_bonds(), rdkit_molecule.GetBonds())
    for bond, rdkit_bond in bonds:
        assert bond.order == rdkit_bond.GetBondTypeAsDouble()
        assert bond.atom1.id == rdkit_bond.GetBeginAtomIdx()
        assert bond.atom2.id == rdkit_bond.GetEndAtomIdx()


def is_equivalent_atom(atom1, atom2):
    return (
        atom1 is not atom2
        and atom1.id == atom2.id
        and atom1.charge == atom2.charge
        and atom1.__class__ is atom2.__class__
    )


def is_equivalent_bond(bond1, bond2):
    return (
        bond1 is not bond2
        and bond1.__class__ is bond2.__class__
        and bond1.order == bond2.order
        and bond1.periodicity == bond2.periodicity
        and is_equivalent_atom(bond1.atom1, bond2.atom1)
        and is_equivalent_atom(bond1.atom2, bond2.atom2)
    )
