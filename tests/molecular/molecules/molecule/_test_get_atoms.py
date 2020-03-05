import stk
import itertools as it
import rdkit.Chem.AllChem as rdkit


def test_get_atoms(tmp_path, case_data):
    # Visual inspection of the molecule can be useful.
    case_data.molecule.write(tmp_path / 'molecule.mol')
    _test_get_atoms(case_data.molecule, case_data.smiles)


def _test_get_atoms(molecule, smiles):
    # The basic idea is that you get a molecule and the canonical
    # smiles for that molecule. You then create a clone "molecule"
    # making sure that atoms have canonical ordering. This conversion
    # will only work correctly when get_atoms() and get_bonds() work
    # correctly - see get_rdkit_mol(). Finally, you compare that atoms
    # of the canonically ordered molecule to an rdkit molecule which
    # also has canonical ordering, and all the atoms should match.

    print(rdkit.MolToSmiles(molecule.to_rdkit_mol()))

    molecule = with_canonical_atom_ordering(molecule)
    expected = rdkit.AddHs(rdkit.MolFromSmiles(smiles))
    for atom1, atom2 in it.zip_longest(
        molecule.get_atoms(),
        expected.GetAtoms(),
    ):
        is_equivalent_atom(atom1, atom2)


def with_canonical_atom_ordering(molecule):
    return stk.BuildingBlock(
        smiles=rdkit.MolToSmiles(get_rdkit_mol(molecule))
    )


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
            order=rdkit.BondType(bond.get_order())
        )
    return rdkit_mol.GetMol()


def is_equivalent_atom(atom, rdkit_atom):
    assert atom.get_atomic_number() == rdkit_atom.GetAtomicNum()
    assert atom.get_id() == rdkit_atom.GetIdx()
    assert atom.get_charge() == rdkit_atom.GetFormalCharge()
